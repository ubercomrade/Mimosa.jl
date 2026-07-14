# Security

## Safe parsing

Mimosa.jl is designed to safely parse untrusted user input:

- **Size limits**: All parsers enforce maximum file sizes and matrix dimensions
- **Strict validation**: Invalid data produces typed errors, not silent
  corruption
- **No `eval`**: No code execution during parsing
- **No unsafe deserialization**: No `pickle`, `joblib`, or Julia
  `Serialization` in the core library
- **No global state**: `using Mimosa` does not touch the filesystem, create
  directories, launch threads, or modify global settings

Portable model and null bundles additionally limit manifests, blob sizes,
array counts, ranks, dimensions, element counts and total declared allocation.
Manifest paths are relative and resolved targets must stay below the bundle
root; checksum and NPY metadata validation happens before array allocation.

## Error handling

```julia
try
    model = readmodel("user_input.meme")
catch e
    if e isa ModelFormatError
        # User provided an invalid file
        println("Invalid model file: $(e.path) — $(e.message)")
    elseif e isa ModelDimensionError
        # Matrix has wrong dimensions
        println("Dimension error: $(e.message)")
    end
end
```

## File system safety

- **Cache**: Directory not created on `Cache()` construction; only on first
  `cache_set`. `using Mimosa` never creates directories.
- **Atomic writes**: Model and null storage commit a complete sibling staging
  directory with rename; orphan staging directories are ignored by readers
- **Checksum validation**: Model bundles and cache entries validate SHA-256
  on every access. Corrupted files result in cache misses or typed errors,
  not silent corruption.

## Legacy format handling

Legacy Python `pickle`/`joblib` files may contain arbitrary Python objects.
Mimosa.jl **never** reads these files directly. Conversion is done via
separate Python scripts with an explicit `--trusted-input` flag:

```bash
python scripts/convert_legacy_model.py --trusted-input old_model.pkl --output new_bundle
```

The `--trusted-input` flag must be provided explicitly, serving as a security
guard against accidental deserialization of untrusted data.

## XML safety

The minimal XML parser for Dimont/Slim models:
- Does not resolve external entities
- Does not make network requests
- Does not execute arbitrary code
- Enforces size limits on input files
