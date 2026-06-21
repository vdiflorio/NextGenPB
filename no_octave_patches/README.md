# NextGenPB senza Octave — patch

Patch per compilare/eseguire NextGenPB **senza Octave**, con bimpp costruito
localmente. Dettagli completi, trappole e comandi: `../../BIMPP_NO_OCTAVE_PLAN.md`.

## File

| File | Dove applicarlo | Cosa contiene |
|------|-----------------|---------------|
| `octave_file_io.patch` | dentro `octave_file_io/` | shim dei tipi Octave (`octave_file_io_shim.h`), branch `OFIO_NO_OCTAVE` in `octave_file_io.h`, text-backend in `octave_file_io.cpp`, `build_no_octave.sh` |
| `bimpp.patch` | dentro `bimpp/` | fix `configure.ac`: `AC_SUBST` vuoti nel ramo "Octave assente" (additiva, zero regressioni) → niente stub mkoctfile |
| `bimpp.build_configuration.sh` | copia in `bimpp/build/` | invocazione `configure` no-Octave (prefix `install_local`, `--with-octave-home=/nonexistent`, `-DOFIO_NO_OCTAVE`). `build/` è gitignored in bimpp, quindi NON è nel `git diff` |

## Applicazione

```sh
# 1) octave_file_io
cd octave_file_io && git apply /path/to/octave_file_io.patch

# 2) bimpp
cd bimpp && git apply /path/to/bimpp.patch
cp /path/to/bimpp.build_configuration.sh build/build_configuration.sh
```

## Build

Macro chiave: `-DOFIO_NO_OCTAVE` (seleziona lo shim al posto di `<octave/oct.h>`).
Senza la macro il comportamento resta quello originale (con Octave).

1. `octave_file_io`: `sh build_no_octave.sh` → `install_no_octave/` (lib+header).
2. `bimpp`: `./autogen.sh` (rigenera `configure` dopo la patch), poi
   `build/build_configuration.sh`, `make -C src`, `make -C src install`, e copia
   `bim_config.h`/`bim_config_pre.h` generati in `install_local/include/`.
3. `NextGenPB`: `local_settings.mk` punta a `octave_file_io/install_no_octave` e
   `bimpp/install_local`; `make`.

### Note specifiche del path con `&` (questa macchina)

`…/new_ns&ngpb`: autotools rifiuta il `&` nella working dir → bimpp si builda via
symlink senza `&` e poi si riscrivono gli `install_name` sul path reale con
`install_name_tool` (dyld gestisce il `&` a runtime). Su un path normale niente di
tutto ciò. Sequenza completa in `../../BIMPP_NO_OCTAVE_PLAN.md`.
