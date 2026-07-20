# План оптимизации profile comparison и построения null

## Цели

- Использовать один процесс Julia и параллелизм только через `Threads`.
- Сделать empirical normalization общей для обычного `compare`, `prepare_profile` и `build_null`.
- Уменьшить пиковое потребление памяти, количество аллокаций и нагрузку на GC.
- Явно разделить параллелизм уровня моделей/пар и уровня последовательностей.
- Сохранить статистическое определение null, seed schedule и численные результаты.

## Вне области изменений

- Не уменьшать background, количество foreground-последовательностей или `n_samples`.
- Не использовать pool переиспользуемых shuffled-профилей.
- Не менять `shuffle` и алгоритм выборки пар.
- Не добавлять специализированный PWM scan/kernel.
- Не добавлять process workers, `Distributed`, `pmap` или внешний multiprocessing.

## Текущее состояние

- В Mimosa.jl уже нет process workers: `ThreadedExecution` использует `Threads.@spawn` и общую память процесса.
- CLI-параметр `--jobs` является устаревшим псевдонимом `--threads`, а не числом Julia workers.
- `build_null` распараллеливает пары, но передаёт `SerialExecution()` во внутренние вызовы `prepare_profile`.
- Обычный scalar `compare` использует `execution` для сканирования последовательностей, а vector `compare` — для target-моделей.
- Empirical normalization реализована несколькими путями: `_fit_transform_empirical`, `fit(EmpiricalLogTail(), flatten_bundle(...))` и `normalize_bundle`.
- Для отдельного background создаются `bg_raw`, объединённая копия scores, дополнительная `Float32`-копия в `fit`, `sortperm`, `unique_scores` и `log_tail`.

## Этап 1. Зафиксировать эталонную корректность и производительность

До изменения алгоритма добавить детерминированные тестовые эталоны для:

- `LogTailTable.scores` и `LogTailTable.log_tail`;
- нормализованных forward/reverse tracks;
- anchors;
- обычного motif `compare`;
- `build_null(shuffle=true)` — `raw_scores`, `pairs` и sampling metadata;
- serial и threaded режимов.

Основное требование — побитовое равенство старого и нового результата там, где вход представлен `Float32`:

```julia
new_table.scores == reference_table.scores
new_table.log_tail == reference_table.log_tail
new_result.raw_scores == reference_result.raw_scores
```

Добавить baseline-бенчмарки:

1. `ScoreProfile` с self-normalization.
2. Motif profile при `background === sequences`.
3. Motif profile с отдельным большим background.
4. Небольшой `build_null(shuffle=true)`.
5. Serial, outer-threaded и scan-threaded варианты.

Измерять elapsed time, allocations, allocated bytes, peak RSS и comparisons/hour. Peak RSS измерять в отдельном Julia-процессе, чтобы GC и предыдущие бенчмарки не искажали результат.

Затрагиваемые файлы:

- `test/unit/test_profiles.jl`;
- `test/unit/test_null_distribution.jl`;
- `test/unit/test_parallel.jl`;
- `benchmark/runbenchmarks.jl`.

## Этап 2. Ввести единый pipeline empirical normalization

Оставить одну каноническую внутреннюю точку входа, например:

```julia
_fit_normalize_empirical(
    raw::StrandPair;
    calibration::StrandPair=raw,
    execution::ExecutionPolicy=SerialExecution(),
)
```

Контракт функции:

1. Построить `LogTailTable` строго из `calibration`.
2. Применить таблицу к `raw`.
3. Вернуть `(table, normalized_bundle)`.
4. Одинаково обработать self-normalization и отдельный background.
5. Корректно обработать как независимые цепи, так и случай `forward === reverse`.

Все рабочие пути должны использовать только этот pipeline:

```text
compare
  -> prepare_profile
      -> _fit_normalize_empirical

build_null
  -> prepare_profile
      -> _fit_normalize_empirical
```

Неиспользуемый `_resolve_profile_bundle` удалить, чтобы не осталось второго независимого алгоритма нормализации.

Основные файлы:

- `src/profiles/normalization.jl`;
- `src/profiles/alignment.jl`;
- `src/comparison/profile_comparison.jl`;
- `src/statistics/null_distribution.jl`.

## Этап 3. Заменить `sortperm` на consuming in-place fit

Добавить внутреннюю функцию:

```julia
_fit_empirical_table!(workspace::Vector{Float32}) -> LogTailTable
```

Она получает принадлежащий pipeline рабочий массив и имеет право его изменить:

1. Выполнить `sort!(workspace; rev=true)`.
2. Последовательно найти группы одинаковых `Float32`.
3. Записать unique scores в начало того же массива.
4. Рассчитать прежнее значение `-log10(cumulative_count / total_count)`.
5. Выполнить `resize!(workspace, n_unique)`.
6. Вернуть workspace как `LogTailTable.scores` и отдельно выделенный `log_tail`.

Это должно удалить:

- `sortperm::Vector{Int}`, пропорциональный размеру background;
- повторную `Float32`-копию внутри `fit`;
- отдельную аллокацию `unique_scores`.

Публичный `fit` не должен мутировать пользовательский вход:

```julia
function fit(::EmpiricalLogTail, scores::AbstractVector{<:Real})
    workspace = Float32.(scores)
    return _fit_empirical_table!(workspace)
end
```

Внутренний pipeline должен передавать уже принадлежащий ему workspace напрямую.

### Требование численной совместимости

Алгоритм остаётся статистически и численно прежним:

- сортировка score по убыванию;
- точное сравнение значений `Float32`;
- прежний cumulative tail count;
- прежний `Float64`-расчёт отношения и логарифма с преобразованием результата в `Float32`;
- прежнее правило lookup между соседними score.

`normalization_version = "empirical-log-tail-v1"` и версию cache менять только если тесты обнаружат изменение результата. При отсутствии побитового равенства изменение не принимать без отдельного анализа совместимости.

## Этап 4. Оптимизировать calibration workspace и время жизни памяти

Добавить внутреннюю функцию, которая сразу копирует bundle в один массив окончательного размера:

```julia
_empirical_workspace(bundle::StrandPair) -> Vector{Float32}
```

Требования:

- одна итоговая аллокация `Vector{Float32}`;
- заполнение через `copyto!`;
- отсутствие цепочки `vcat` и последующей конвертации;
- сохранение текущей семантики для `forward === reverse`;
- проверка переполнения при расчёте общей длины.

Для отдельного background организовать время жизни объектов так:

```text
foreground raw
background raw
calibration workspace
fit table
освобождение background raw
normalization foreground
anchors
```

Подготовку таблицы из background вынести в небольшую отдельную функцию, чтобы компилятор мог раньше завершить время жизни `bg_raw`. После реализации проверить это по peak RSS, а не полагаться только на оценку аллокаций Julia.

## Этап 5. Потоковое применение готовой таблицы

Сам fit таблицы оставить последовательным: стандартный `sort!` не является многопоточным, а новая зависимость для parallel sort в этот этап не входит.

Применение готовой таблицы к foreground можно распараллелить:

```julia
normalize_bundle(table, bundle; execution=SerialExecution())
```

Делить массивы на крупные непрерывные диапазоны. Поэлементное преобразование независимо, поэтому такой threading не меняет результат и не требует синхронизации.

В режиме outer threading применение таблицы внутри каждой пары должно оставаться serial, чтобы не создавать вложенный параллелизм и oversubscription.

## Этап 6. Разделить уровни параллелизма в `build_null`

Расширить API двумя явными политиками:

```julia
build_null(
    models;
    execution=ThreadedExecution(2),
    scan_execution=SerialExecution(),
    ...,
)
```

Семантика:

- `execution` — верхний уровень, то есть пары null;
- `scan_execution` — последовательности внутри одного `prepare_profile`, включая применение normalization table.

Поддержать два режима.

### Параллельные пары

```julia
execution=ThreadedExecution(2)
scan_execution=SerialExecution()
```

Этот режим допускает несколько одновременно готовящихся пар и подходит при достаточном объёме RAM.

### Параллельное сканирование одной пары

```julia
execution=SerialExecution()
scan_execution=ThreadedExecution(6)
```

Этот режим ограничивает число одновременно живущих больших профилей, но использует несколько потоков при сканировании последовательностей и применении нормализации.

### Запрет вложенного threading

Если обе политики являются `ThreadedExecution` с фактическим числом задач больше одного, выбрасывать понятный `ArgumentError`. Не следует молча превращать внутренний уровень в serial только через `_in_parallel_region()`.

Внутренняя защита `_in_parallel_region()` остаётся последней гарантией от случайного вложенного параллелизма.

## Этап 7. Унифицировать уровни параллелизма обычного `compare`

Применить те же понятия к profile comparisons:

- `execution` — верхний доступный уровень: target-модели для vector comparison;
- `scan_execution` — последовательности внутри `prepare_profile`.

Для scalar comparison верхнего уровня коллекции нет, поэтому рекомендуемый режим:

```julia
execution=SerialExecution()
scan_execution=ThreadedExecution(n)
```

Для vector comparison:

```julia
execution=ThreadedExecution(n)
scan_execution=SerialExecution()
```

Обратную совместимость сохранить:

- существующие вызовы с одним `execution` продолжают работать;
- текущие defaults не меняют результаты и уровень параллелизма;
- новый `scan_execution` по умолчанию равен `SerialExecution()` там, где существует outer level;
- изменение поведения scalar overload допускается только с явной документацией и тестами.

Главное инвариантное требование: `compare` и `build_null` используют один `prepare_profile` и один empirical normalization pipeline.

## Этап 8. Оставить только threads в API и документации

Process workers в библиотеку не добавлять. Терминологию привести к одному варианту:

- основной CLI-флаг — `--threads`;
- `--jobs` либо удалить в следующем breaking release, либо временно принимать с deprecation warning;
- явно документировать один процесс Julia и общую память;
- показать запуск runtime с нужным количеством потоков:

```bash
julia --threads=6 -m Mimosa build-null ... --threads=6
```

Внешний рекомендуемый workflow также должен использовать один Julia-процесс. `addprocs`, `pmap` и параллельный запуск нескольких экземпляров не входят в поддерживаемую схему выполнения null.

Затрагиваемые файлы:

- `src/cli.jl`;
- `docs/src/cli.md`;
- `docs/src/api.md`;
- `docs/src/quickstart.md`.

## Этап 9. Тестирование параллелизма и type stability

Добавить следующую матрицу тестов:

| Операция | Outer | Scan | Ожидаемый результат |
|---|---|---|---|
| scalar `compare` | serial | serial | reference |
| scalar `compare` | serial | threaded | идентично reference |
| vector `compare` | threaded | serial | идентично reference |
| `build_null` | serial | serial | reference |
| `build_null` | threaded | serial | идентично reference |
| `build_null` | serial | threaded | идентично reference |
| любая операция | threaded | threaded | `ArgumentError` |

Дополнительно проверить:

- каждый work item обрабатывается ровно один раз;
- порядок `raw_scores` и `pairs` не меняется;
- seed schedule не меняется;
- исключения из threaded scan корректно пробрасываются;
- число активных задач ограничено политикой и `Threads.nthreads()`;
- public `fit` не мутирует вход;
- consuming fit используется только с внутренним принадлежащим ему workspace;
- empty input, single value, duplicates, signed zero и независимые strand bundles;
- `@code_warntype` для нового normalization pipeline, `prepare_profile` и обоих режимов `build_null` не показывает новых динамически типизированных hot paths.

## Этап 10. Финальная проверка производительности

Повторить baseline на одинаковых данных и runtime-настройках.

Обязательные сравнения:

- старая и новая нормализация с отдельным background;
- serial до/после;
- outer threads: 1, 2, 3, 4;
- scan threads: 1, 2, 3, 4, 6;
- peak RSS для representative promoter background;
- allocations на одну подготовку shuffled-профиля;
- comparisons/hour для representative `build_null`.

Основной критерий выбора режима — throughput при приемлемом peak RSS, а не средняя загрузка CPU.

## Критерии готовности

Модификация считается завершённой, если:

- в исходниках нет `Distributed` и process workers;
- обычное сравнение и null используют один normalization pipeline;
- отдельный background fit не создаёт `sortperm`, пропорциональный числу scores;
- отсутствуют лишняя конвертирующая копия и отдельная аллокация `unique_scores`;
- public `fit` не мутирует пользовательский вход;
- serial, outer-threaded и scan-threaded результаты идентичны;
- seed schedule, shuffle и null sampling не изменены;
- не добавлен специализированный PWM-код;
- не изменены размеры статистических выборок;
- одновременный outer и inner threading явно запрещён;
- peak RSS на representative background снизился;
- serial performance не ухудшился;
- новый код type-stable и покрыт тестами;
- документация однозначно описывает threads и выбор одного уровня параллелизма.

## Порядок реализации

1. Эталонные correctness-тесты и baseline-бенчмарки.
2. Единый empirical normalization pipeline.
3. Consuming in-place fit без `sortperm`.
4. Calibration workspace без промежуточных копий.
5. Потоковое применение готовой таблицы.
6. Разделение outer и scan execution в `build_null`.
7. Унификация execution API обычного `compare`.
8. Тесты type stability и полной матрицы параллелизма.
9. CLI и документация только в терминах threads.
10. Финальные бенчмарки времени, аллокаций и peak RSS.
