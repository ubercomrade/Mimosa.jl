# План развития API расширения моделей Mimosa.jl

Статус: проект плана, реализация не начата.

Дата ревизии: 2026-07-14.

## 1. Цель

Сделать пользовательские модели полноценными участниками библиотечных
workflow Mimosa.jl без изменения исходников пакета. Для базового сравнения
пользователь должен определить тип модели, ее имя, длину мотива и одно ядро,
вычисляющее оценки обеих цепей. Геометрия окна, пакетное сканирование,
параллелизм, нормализация, подготовка профиля и сравнение должны предоставляться
Mimosa.jl.

Целевой путь использования:

```julia
import Mimosa

struct MyModel{P} <: Mimosa.AbstractMotifModel
    name::String
    parameters::P
    length::Int
end

Mimosa.modelname(model::MyModel) = model.name
Mimosa.motif_length(model::MyModel) = model.length

function Mimosa.scan_kernel!(
    forward::AbstractVector{Float32},
    reverse::AbstractVector{Float32},
    model::MyModel,
    sequence::AbstractVector{UInt8},
    n_positions::Int,
)
    # Вход уже проверен безопасной публичной границей Mimosa.scan/scan!.
    # Заполнить forward[1:n_positions] и reverse[1:n_positions].
    return (forward, reverse)
end

sequences = Mimosa.readsequences("sequences.fasta")
result = Mimosa.compare(MyModel(...), built_in_model, sequences)
```

Для модели с левым или правым контекстом дополнительно определяются только
ненулевые границы:

```julia
Mimosa.left_context(model::MyModel) = 2
Mimosa.right_context(model::MyModel) = 1
```

## 2. Не цели

- Не вводить глобальный изменяемый registry моделей или форматов.
- Не загружать Julia-код, типы или функции из model bundle.
- Не использовать `eval`, Julia `Serialization` или иные небезопасные способы
  восстановления пользовательского типа.
- Не делать внутренние traits оптимизированных higher-order ядер обязательной
  частью внешнего API.
- Не менять Float32-порядок вычислений, tie-breaking, координатные соглашения
  или существующие результаты встроенных моделей.
- Не обещать автоматическое обнаружение стороннего формата стандартным CLI:
  сторонняя интеграция формата сначала является library-only API.
- Не менять model/null/cache format versions, пока фактически не изменена их
  схема. Любое такое изменение требует отдельного решения о миграции.

## 3. Геометрический контракт

### 3.1. Публичные определения

Геометрия модели задается следующими функциями:

```julia
motif_length(model)       # длина возвращаемого сайта
left_context(model) = 0   # число баз слева от сайта, нужных для одной оценки
right_context(model) = 0  # число баз справа от сайта, нужных для одной оценки
```

Mimosa.jl вычисляет производные величины:

```julia
window_size(model) =
    left_context(model) + motif_length(model) + right_context(model)

npositions(model, sequence_length) =
    max(sequence_length - window_size(model) + 1, 0)

site_start_offset(model) = left_context(model)
```

`motif_length` всегда означает число баз в возвращаемом сайте, а не число
столбцов произвольного внутреннего представления. Scan position обозначает
начало полного окна. Начало сайта в однобазной системе Julia равно
`scan_position + left_context(model)`.

`left_context` и `right_context` задаются относительно возрастающих координат
исходной последовательности, а не относительно ориентации мотива. Forward и
reverse score с одним индексом относятся к одному физическому полному окну и к
одному физическому участку сайта внутри него. Поэтому смещение начала сайта
равно `left_context` для обеих цепей; reverse kernel отвечает за правильную
ориентацию вычисления оценки. Если будущей модели нужны разные физические
участки сайта для двух ориентаций, это не покрывается данным контрактом и
требует отдельного ADR, а не неявной перестановки left/right context.

### 3.2. Отображение встроенных моделей

| Модель | `motif_length` | `left_context` | `right_context` |
|---|---:|---:|---:|
| PWM | `length(model)` | `0` | `0` |
| SiteGA | `model.motif_length` | `0` | `0` |
| BaMM | `model.motif_length` | `model.order` | `0` |
| Dimont | `model.motif_length` | `model.span` | `0` |
| Slim | `model.motif_length` | `model.span` | `0` |

`order` и `span` остаются терминами конкретных типов и не входят в общий
контракт `AbstractMotifModel`. Внутренний `context_length` перестает быть
внешней точкой расширения. На время миграции он может делегировать в
`left_context`, если это необходимо для существующих higher-order ядер.

### 3.3. Валидация геометрии

- `motif_length(model)` должен возвращать положительный `Integer`, приводимый к
  `Int` без потери значения.
- `left_context(model)` и `right_context(model)` должны возвращать
  неотрицательные `Integer`.
- Сложение границ и длины должно проверяться на переполнение до выделения
  памяти.
- `npositions` должен возвращать ноль для слишком короткой и пустой
  последовательности.
- Все строки batch должны использовать одну геометрию модели и сохранять
  исходный порядок, включая пустые строки результата.
- Извлекаемый диапазон сайта должен целиком лежать внутри исходной
  последовательности для обеих ориентаций.
- Reverse-complement координаты должны сохранить существующее однобазное
  включительное соглашение библиотеки.
- Для асимметричного контекста forward и reverse score одного индекса должны
  ссылаться на одинаковый физический участок сайта; reverse extraction меняет
  только ориентацию возвращаемых баз.

## 4. Минимальный контракт модели

### 4.1. Обязательные методы для сравнения

Для capability `:compare` обязательны только:

```julia
modelname(model::MyModel)::AbstractString
motif_length(model::MyModel)::Integer
scan_kernel!(forward, reverse, model::MyModel, sequence, n_positions)
```

`left_context` и `right_context` имеют общий default `0`. Пользователь
переопределяет их только для модели, действительно использующей контекст.

`modelname` должен возвращать непустое стабильное имя экземпляра. Алгоритмы не
должны требовать наличия поля `name`.

`scan_kernel!` является стабильной точкой расширения, а не самостоятельной
безопасной границей для произвольного ввода. Перед его вызовом Mimosa.jl должна:

- проверить DNA-коды и геометрию последовательности;
- вычислить и проверить `n_positions`;
- проверить тип, длины и отсутствие aliasing выходных буферов;
- гарантировать, что ядру доступны ровно `n_positions` элементов результата;
- не разрешать ядру менять модель или последовательность.

Ядро заполняет обе цепи и возвращает `(forward, reverse)`. Канонический тип
сканирующих значений остается `Float32`.

### 4.2. Производные strand-операции

Mimosa.jl предоставляет корректные fallback-реализации:

- `scan_both!` вызывает `scan_kernel!` после валидации;
- `scan_forward!` использует forward-часть парного ядра;
- `scan_reverse!` использует reverse-часть парного ядра;
- `best_hits!` вычисляет поэлементный максимум в установленном порядке;
- allocating `scan` выделяет правильные плоские или ragged outputs;
- batch `scan` использует существующий serial/threaded scheduler.

Fallback может вычислять обе цепи даже для одностороннего запроса. Модель может
опционально переопределить безопасно задокументированные специализированные
kernels после benchmark. Эти методы являются ускорениями, а не условием
совместимости.

### 4.3. Встроенные оптимизированные модели

PWM, BaMM, SiteGA, Dimont и Slim сохраняют текущие специализированные ядра и
численный порядок. Общий fallback не должен насильно заменять быстрые пути.
Внутренние `kmer`, rolling-code, `scan_width`, matrix layout и сходные traits
остаются деталями Mimosa.jl.

`is_scannable` следует исключить из обязательного контракта. Само наследование
от `AbstractMotifModel` означает поддержку сканирования. На переходный период
функцию можно сохранить как deprecated compatibility shim, возвращающий
результат проверки интерфейса.

### 4.4. Необязательные capabilities

- `scorebounds(model)` нужен для inspect/экспорта форматов, которым необходимы
  границы, но не для `compare`.
- `model_fingerprint(model)` нужен для cache и null bundles.
- Typed reader/writer нужен только для файлового ввода-вывода.
- Специализированные scan kernels нужны только для производительности.

Отсутствие необязательного capability не должно мешать прямому
`scan`/`prepare_profile`/`compare`.

## 5. Устранение структурного duck typing

### 5.1. Новые публичные accessors

Добавить и экспортировать:

```julia
modelname(source::AbstractProfileSource)
motif_length(model::AbstractMotifModel)
left_context(model::AbstractMotifModel) = 0
right_context(model::AbstractMotifModel) = 0
window_size(model::AbstractMotifModel)
site_start_offset(model::AbstractMotifModel)
model_fingerprint(source::AbstractProfileSource)
```

`window_size` и `site_start_offset` получают общие реализации из раздела 3.
Методы встроенных типов должны подтверждать эти формулы тестами. Прямые
переопределения производных функций допустимы только как временная миграция или
для модели, геометрия которой не представима контрактом; второй случай требует
отдельного проектного решения.

### 5.2. Замена доступа к полям

Во всех обобщенных алгоритмах заменить:

- `model.name` на `modelname(model)`;
- `model.motif_length` и `length(model)` на `motif_length(model)`, если речь о
  длине сайта;
- `model.order`/`model.span` на `left_context` или конкретный внутренний dispatch;
- `model.representation` на capability-specific dispatch, если
  код не является реализацией конкретного встроенного типа.

Области обязательного аудита:

- `src/comparison/` и `src/profiles/alignment.jl`;
- `src/statistics/null_distribution.jl`;
- `src/cache/cache.jl`;
- `src/sites/sites.jl`;
- `src/io/model_storage.jl`;
- `src/cli.jl`;
- precompile workload и JSON/result adapters.

Доступ к полям остается нормальным внутри конструктора, ядра или codec
конкретного встроенного типа. Запрещается только предположение, что сторонний
подтип имеет те же поля.

### 5.3. Fingerprint

Текущая generic-функция fingerprint распознает встроенные поля через `isa` и
отклоняет пользовательские модели. Ее необходимо разделить:

- публичный `model_fingerprint(source)` возвращает стабильную SHA-256 строку;
- встроенные типы сохраняют текущие байтовые представления и cache keys;
- сторонняя модель реализует метод только при запросе cache/null capability;
- сравнение без cache не вызывает fingerprint;
- fingerprint обязан включать тип/версию научной модели и все параметры,
  влияющие на оценки;
- имя экземпляра включается только если это соответствует уже закрепленному
  контракту cache; менять существующие ключи без миграции нельзя.

Не следует предоставлять generic fingerprint через отражение полей: порядок и
представление полей не являются стабильным форматом.

## 6. Профили как простой адаптер

### 6.1. Целевые конструкторы

Предоставить ergonomic-конструкторы:

```julia
ScoreProfile(name, rows)
ScoreProfile(name, forward_rows, reverse_rows)
ScoreProfile(name, forward::RaggedArray, reverse::RaggedArray)
```

Первый конструктор создает симметричный профиль и сохраняет текущее поведение.
Второй и третий сохраняют независимые forward/reverse tracks. Это рекомендуемый
путь для внешнего сканера, которому не нужны `AbstractMotifModel`, site
extraction и встроенное пакетное сканирование.

### 6.2. Представление и совместимость

Сделать `ScoreProfile` параметрическим по конкретным типам двух `RaggedArray`,
чтобы избежать `Union` и абстрактных полей. Текущее поле `scores` можно оставить
как forward track на один deprecation cycle либо заменить публичными accessors:

```julia
forward_scores(profile)
reverse_scores(profile)
profile_bundle(profile)
```

Двухаргументный `ScoreProfile(name, scores)` должен продолжать работать и
использовать один и тот же immutable-by-contract ragged объект для обеих цепей,
если это безопасно.

### 6.3. Валидация профилей

- Имя должно быть непустым.
- Значения приводятся к `Float32` на границе конструктора либо отклоняются по
  единому задокументированному правилу.
- Все значения должны быть конечными.
- Forward и reverse должны иметь одинаковое число строк.
- Соответствующие forward/reverse строки должны иметь одинаковую длину.
- Пустые строки и исходный порядок должны сохраняться.
- Конструктор из строк может использовать временный `Vector{Vector}` как
  удобную холодную границу, но хранение остается плоским `RaggedArray`.
- `prepare_profile` и `compare` не должны считать обе цепи одинаковыми для
  strand-aware профиля.

### 6.4. Чтение score-файлов

Текущий формат `read_scores` остается симметричным для совместимости. Для
strand-aware файла нужно либо задокументировать отдельный bounded format, либо
предоставить два явных файла/потока. Изменение синтаксиса существующего файла
требует format version и негативных security tests.

## 7. Обобщение связанных workflows

### 7.1. Сравнение и подготовка

Существующие методы `prepare_profile(::AbstractMotifModel, ...)` и
`compare(::AbstractMotifModel, ...)` должны зависеть только от минимального
контракта и `modelname`. Проверить следующие комбинации:

- custom/custom;
- custom/built-in и built-in/custom;
- prepared custom/raw target;
- raw custom/prepared target;
- one-to-many с custom query и смешанным набором targets;
- serial и explicit threaded execution.

Mixed `ScoreProfile`/motif путь может оставаться запрещенным до явного решения
о нормализации, но сообщение об ошибке должно предлагать подготовить обе стороны
как profiles.

### 7.2. Site selection и PFM reconstruction

Заменить дублированные публичные dispatch для `PWM` и
`AbstractHigherOrderMotif` общими методами для `AbstractMotifModel`.

Общий `selectsites` должен использовать только:

- `scan(...; strands=BothStrands())`;
- `motif_length`;
- `left_context`/`site_start_offset`;
- общую логику selector и сортировки.

Общий `reconstruct_pfm` должен извлекать участок длины `motif_length`, начиная
после `left_context`, и не включать ни левый, ни правый контекст. Проверить обе
ориентации и модели с одновременно ненулевыми left/right context.

### 7.3. Null distributions и cache

`build_null` должен использовать `modelname` и минимальный compare capability.
Fingerprint требуется только при создании артефакта, где он входит в
совместимость или cache key. Ошибка отсутствующего fingerprint должна называть
конкретный capability и метод, который нужно реализовать.

Не расширять типы сохраненного null bundle пользовательскими Julia type names.
Схема null version 2 остается profile-only.

## 8. Отделение файловых форматов от моделей

### 8.1. Typed format dispatch

Ввести неизменяемые типы форматов:

```julia
abstract type AbstractModelFormat end

struct AutoModelFormat <: AbstractModelFormat end
struct MemeFormat <: AbstractModelFormat end
struct PFMFormat <: AbstractModelFormat end
struct BaMMFormat <: AbstractModelFormat end
struct SiteGAFormat <: AbstractModelFormat end
struct DimontFormat <: AbstractModelFormat end
struct SlimFormat <: AbstractModelFormat end
struct PortableModelBundle <: AbstractModelFormat end
```

Основные точки dispatch:

```julia
readmodel(path, format::AbstractModelFormat; kwargs...)
writemodel(path, model, format::AbstractModelFormat; kwargs...)
```

Существующий keyword `format=:auto` остается совместимой I/O-границей и
преобразуется во встроенный typed format. Строки и symbols не попадают в
научные kernels.

Сторонний пакет определяет собственный тип формата и более специфичные методы,
не переопределяя широкий `readmodel(path::AbstractString; ...)` и не меняя
исходники Mimosa.jl.

### 8.2. Auto-detection

Встроенный `AutoModelFormat` обнаруживает только форматы, известные Mimosa.jl.
Не вводить глобальную регистрацию расширений при `using Mimosa`. Сторонний
пользователь передает формат явно:

```julia
model = readmodel(path, MyModelFormat())
```

Внешний пакет или приложение может реализовать собственный auto-detection за
пределами Mimosa.jl.

### 8.3. Portable bundles

Текущий portable model bundle продолжает принимать только встроенные kinds,
пока не принят отдельный versioned codec protocol. Нельзя сохранять
пользовательский Julia type name и затем динамически его конструировать.

Если generic custom bundle понадобится, отдельный проект должен определить:

- глобально уникальный codec/provider identifier;
- versioned schema параметров;
- bounded declarations и лимиты до allocation;
- checksum и path/symlink validation;
- явный decoder, уже доступный в процессе;
- поведение при отсутствии decoder;
- миграции и влияние на model/cache format versions.

### 8.4. CLI

Стандартный CLI не должен принимать имя Julia-модуля или выражение для
динамической загрузки. Безопасные варианты будущей интеграции:

- отдельный CLI, предоставляемый пакетом модели;
- предварительное преобразование в поддерживаемый portable формат;
- strand-aware score profiles как нейтральный обменный формат.

## 9. Проверяемый внешний контракт

### 9.1. Runtime validation

Добавить и экспортировать:

```julia
validate_model(model; capability=:compare)
```

Symbols допустимы на этой пользовательской границе. Внутри capability следует
преобразовать в небольшой concrete type либо выполнить только холодную
валидацию, не включая string dispatch в hot path.

Поддерживаемые capabilities:

| Capability | Проверки |
|---|---|
| `:compare` | `AbstractMotifModel`, имя, геометрия, pair kernel |
| `:sites` | `:compare` плюс корректное смещение и диапазон сайта |
| `:cache` | `:compare` плюс стабильный fingerprint |
| `:storage_read` | typed reader для заданного format |
| `:storage_write` | typed writer для пары model/format |

Функция должна возвращать модель либо `nothing` при успехе и выбрасывать
отдельный `ModelInterfaceError <: MimosaError` при нарушении. Ошибка должна
содержать capability, тип модели, отсутствующие методы и обнаруженные неверные
значения. Не следует полагаться на перехват произвольного `MethodError` из
глубины вычисления.

### 9.2. Граница валидации

Полную interface validation выполнять один раз в публичных `scan`,
`prepare_profile`, `compare`, `selectsites`, cache/storage entry points.
Внутренние циклы и worker tasks получают уже проверенную модель и не повторяют
динамическую проверку на каждой позиции.

Worker exception из пользовательского ядра должен распространяться вызывающей
задаче. Частично заполненный результат не должен возвращаться.

### 9.3. Downstream contract test

В `test/downstream/` определить модель в отдельном модуле, импортирующем Mimosa
как обычная зависимость. Тест не должен обращаться к `Mimosa._private_name` и не
должен иметь поля `representation`, `weights`, `order` или `span`.

Обязательные проверки:

- минимальная модель без контекста;
- модель с ненулевыми left и right context;
- `validate_model` success и понятные failures для каждого метода;
- single scan для четырех strand policies;
- empty, too-short, exact-window и longer sequences;
- ragged batch с пустыми строками;
- serial/threaded exact equivalence и сохранение порядка;
- worker exception propagation;
- nested threaded fallback;
- scalar, prepared и one-to-many comparison;
- custom/built-in comparison в обоих порядках;
- site coordinates и PFM extraction без контекста;
- profile adapter с разными forward/reverse tracks;
- отсутствие требования fingerprint для обычного compare;
- явная ошибка cache/null при отсутствии fingerprint;
- отсутствие файловой активности при `using Mimosa`.

## 10. Этапы реализации

### Этап 0. Зафиксировать контракт

- Добавить ADR с определениями scan position, motif site, left/right context и
  reverse coordinates.
- Обновить `docs/src/data_layout.md` и `numerical_compatibility.md` до изменения
  dispatch.
- Добавить characterization tests текущих встроенных моделей: длины строк,
  значения, offsets и ориентации.
- Подтвердить, что переход на производные geometry functions не меняет
  существующие результаты.

Критерий завершения: формулы раздела 3 закреплены тестами для всех пяти
встроенных типов.

### Этап 1. Ввести accessors и новую геометрию

- Добавить `modelname`, `left_context`, `right_context` и общие производные
  методы.
- Реализовать accessors для встроенных моделей и `ScoreProfile`.
- Перевести built-in geometry tests на публичные функции.
- Сохранить временные compatibility shims для `context_length`, прямых
  `window_size` и `site_start_offset`, где это нужно внутренним ядрам.
- Экспортировать только стабильные внешние функции; внутренние traits не
  экспортировать.

Критерий завершения: built-in serial/threaded scan и site tests дают прежние
точные результаты.

### Этап 2. Ввести минимальное scan extension point

- Добавить `scan_kernel!` и безопасную wrapper-границу.
- Вывести четыре strand policies из парного ядра.
- Сохранить specialized built-in dispatch.
- Добавить `ModelInterfaceError` и `validate_model(:compare)`.
- Добавить первую внешнюю модель в downstream tests.
- Проверить type stability built-in hot paths и отсутствие новых allocations в
  существующих specialized kernels.

Критерий завершения: внешний тип с тремя обязательными методами проходит
single/batch/prepared/compare tests без доступа к private API.

### Этап 3. Удалить структурный duck typing

- Перевести comparison, profiles, null, cache, sites и CLI на accessors.
- Разделить минимальный compare и optional fingerprint/scorebounds capabilities.
- Сохранить текущие fingerprints встроенных моделей побайтно.
- Добавить негативный тест модели с другим набором полей.

Критерий завершения: поиск по обобщенным алгоритмам не находит прямых обращений
к полям, принадлежащим только встроенным моделям.

### Этап 4. Расширить ScoreProfile

- Изменить представление на два concrete ragged tracks с совместимым
  симметричным конструктором.
- Добавить ergonomic constructors и accessors.
- Обновить normalization, prepare и compare для независимых цепей.
- Добавить форматный план для strand-aware score input либо явно оставить его
  constructor-only в первой версии.
- Проверить, что старые симметричные результаты не изменились.

Критерий завершения: пользователь может сравнить независимые forward/reverse
выходы внешнего сканера без определения `AbstractMotifModel`.

### Этап 5. Обобщить workflows

- Объединить PWM/higher-order site dispatch в методы `AbstractMotifModel`.
- Перевести PFM reconstruction на новую геометрию.
- Обобщить null construction и one-to-many paths.
- Добавить tests с right context и mixed model collections.

Критерий завершения: custom model получает sites и reconstruction без
дополнительных методов помимо геометрии и сканирования.

### Этап 6. Ввести typed format API

- Добавить `AbstractModelFormat` и built-in concrete format types.
- Перенести чтение/запись на dispatch по format object.
- Оставить keyword-symbol compatibility layer.
- Добавить downstream format type и reader/writer tests.
- Не расширять portable bundle kinds в рамках этого этапа.

Критерий завершения: внешний пакет добавляет явный формат без изменения Mimosa,
глобального registry и type piracy.

### Этап 7. Документация и deprecation cleanup

- Полностью переписать `docs/src/extending.md` вокруг минимального контракта.
- Добавить два законченных примера: custom model и external score adapter.
- Обновить API reference, models, storage, CLI, downstream contract и
  architecture pages.
- Добавить `[Unreleased]` записи в changelog для новых exports и deprecations.
- Удалять shims только в заявленном breaking release.

Критерий завершения: документационный пример копируется в отдельный downstream
module и проходит без внутренних имен.

## 11. Совместимость и выпуск

- `window_size`, `npositions` и `site_start_offset` остаются публичными, но
  становятся производными по умолчанию.
- `order` и `span` не удаляются из встроенных structs или bundle manifests.
- Старый `ScoreProfile(name, scores)` сохраняется.
- Символьный `readmodel(...; format=:...)` сохраняется как boundary adapter.
- Existing built-in model/null/cache schemas не меняются только из-за нового
  extension API.
- Deprecated low-level API должен выдавать `Base.depwarn` только на внешнем
  вызове, не в hot loops.
- Новые exports требуют Aqua ambiguity/piracy checks и export documentation.
- Julia 1.10 остается минимальной поддерживаемой версией.

## 12. Проверка каждого этапа

Сначала запускать узкие тесты затронутой подсистемы, затем downstream contract и
доступную широкую suite. Для каждого этапа обязательны:

```bash
julia --project=Mimosa.jl/test -e \
  'using Mimosa, Test; include("Mimosa.jl/test/unit/test_validation.jl")'

julia --project=Mimosa.jl/test/downstream Mimosa.jl/test/downstream/runtests.jl

JULIA_NUM_THREADS=4 julia --project=Mimosa.jl/test -e \
  'using Mimosa, Test; include("Mimosa.jl/test/unit/test_parallel.jl")'

julia --project=Mimosa.jl/test -e \
  'using JuliaFormatter; @assert format("Mimosa.jl/src"; overwrite=false); @assert format("Mimosa.jl/test"; overwrite=false)'
```

Дополнительно запускать соответствующие `test_profiles.jl`, `test_sites.jl`,
`test_cache.jl`, `test_null_distribution.jl`, `test_model_storage.jl`, CLI,
Aqua и JET tests. Полную suite нельзя объявлять зеленой, пока сохраняется
известный конфликт удаленного root oracle corpus с `test/runtests.jl`.

Для этапов 2 и 5 нужен benchmark до/после на одинаковом public path. Он должен
подтвердить отсутствие регрессии встроенных specialized kernels; generic
fallback пользовательской модели оценивается отдельно и не сравнивается с
ручным specialized kernel как эквивалентный performance contract.

## 13. Итоговый Definition of Done

- Для сравнения внешняя модель реализует три обязательных метода; context
  methods нужны только при ненулевом контексте.
- `window_size`, `npositions` и `site_start_offset` согласованно выводятся из
  новой геометрии.
- Ни один обобщенный workflow не требует полей `name`, `representation`,
  `weights`, `order`, `span` или `motif_length`.
- External model работает в scan, prepared comparison, one-to-many, sites и
  reconstruction через публичный API.
- External score producer работает через strand-aware `ScoreProfile` без типа
  модели.
- Cache/storage capabilities явно отделены от compare capability.
- External format добавляется typed dispatch без изменения Mimosa.jl.
- Interface failures диагностируются до запуска worker tasks.
- Built-in discrete outputs и Float32 values сохраняются; изменение numerical
  contract отсутствует либо оформлено отдельным compatibility decision.
- Документация, downstream tests, formatting, Aqua и JET соответствуют новому
  контракту.
