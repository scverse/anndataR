# Development status

## Introduction

This vignette provides an overview of the current development status of
the *[anndataR](https://bioconductor.org/packages/3.24/anndataR)*
package. It provides details on the current implementation of different
features as well as listing known issues.

## Objects

These tables show the status of the implementation of different
`AnnData` back ends.

### `HDF5AnnData`

| Slot | Getter | Getter test | Setter | Setter test |
|:---|:--:|:--:|:--:|:--:|
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L76) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L366) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L86) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L489) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L227) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L416) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L236) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L500) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L271) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L460) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L280) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L532) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L105) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L376) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L115) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L554) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L169) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L396) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L179) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L677) |
| raw |  |  |  |  |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L315) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L324) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L707) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L249) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L438) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L258) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L516) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L293) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L466) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L302) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L543) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L137) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L386) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L147) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L569) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L198) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L406) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L208) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L692) |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L47) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L359) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L57) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L478) |

### `InMemoryAnnData`

| Slot | Getter | Getter test | Setter | Setter test |
|:---|:--:|:--:|:--:|:--:|
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L75) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L55) | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L79) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L248) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L93) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L16) | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L97) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L148) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L127) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L20) | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L130) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L212) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L149) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L153) |  |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L187) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L191) |  |
| raw |  |  |  |  |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L223) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L226) |  |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L110) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L18) | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L114) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L180) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L138) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L22) | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L141) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L230) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L168) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L172) |  |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L205) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L209) |  |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L57) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L11) | [✅](https://github.com/scverse/anndataR/blob/main/R/InMemoryAnnData.R#L61) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-InMemoryAnnData.R#L119) |

### `ReticulateAnnData`

| Slot | Getter | Getter test | Setter | Setter test |
|:---|:--:|:--:|:--:|:--:|
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L57) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L216) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L74) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L386) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L107) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L206) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L110) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L244) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L137) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L226) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L144) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L367) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L181) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L301) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L184) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L332) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L221) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L313) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L224) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L346) |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L263) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L221) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L266) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L264) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L122) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L211) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L125) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L254) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L159) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L228) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L166) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L373) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L201) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L307) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L204) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L339) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L242) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L319) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L245) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L353) |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L36) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L202) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L39) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L238) |

### `ZarrAnnData`

| Slot | Getter | Getter test | Setter | Setter test |
|:---|:--:|:--:|:--:|:--:|
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L77) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L23) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L86) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L150) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L217) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L71) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L221) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L163) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L255) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L115) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L259) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L199) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L107) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L112) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L225) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L165) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L51) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L170) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L282) |
| raw |  |  |  |  |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L281) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L285) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L313) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L236) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L93) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L240) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L181) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L268) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L121) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L272) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L212) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L136) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L141) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L240) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L191) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L61) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L196) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L298) |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L51) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L16) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L56) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L134) |

## Conversion

These tables show the implementation status of conversion between
`AnnData` and other objects.

### `SingleCellExperiment`

| Slot | From | From test | To | To test |
|:---|:--:|:--:|:--:|:--:|
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L320) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L66) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L209) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L64) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L272) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L24) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L239) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L22) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L271) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L18) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L240) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L18) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L333) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L176) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L347) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L202) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L397) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L90) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L249) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L123) |
| raw |  |  |  |  |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L445) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L145) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L257) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L171) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L296) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L45) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L244) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L43) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L295) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L20) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L245) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L16) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L359) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L183) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L348) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L210) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L421) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L118) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L253) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L147) |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L142) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L208) |  |

### `Seurat`

| Slot | From | From test | To | To test |
|:---|:--:|:--:|:--:|:--:|
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L241) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L73) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L241) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L74) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L189) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L35) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L233) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L33) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L188) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L31) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L266) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L29) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L270) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L124) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L565) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L141) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L342) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L182) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L283) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L99) |
| raw |  |  |  |  |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L415) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L97) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L304) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L113) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L216) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L56) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L250) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L55) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L215) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L29) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L267) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L27) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L300) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L131) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L566) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L149) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L378) |  |  |  |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L138) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L240) |  |

## Known issues

This section lists current known issues in
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)*. Only
certain types of issues are listed here, for additional issues see the
[GitHub issue
tracker](https://github.com/scverse/anndataR/issues?q=sort%3Aupdated-desc+is%3Aissue+is%3Aopen).

### Issue: converted sce object has dimnames(), whilst the original anndata does not.

- Affected backend: `to_SCE`
- Affected slot(s): `obsm`, `varm`
- Affected dtype(s): `pca`
- Probable cause: convert
- To investigate: TRUE
- To fix: FALSE

#### Error message

    sampleFactors(reducedDims(sce)$pca) (`actual`) not equal to ad$obsm[["X_pca"]] (`expected`).
    `dimnames(actual)` is a list `dimnames(expected)` is absent

#### Proposed solution

Investigate if this is a problem or not.

## Session info

``` r

sessionInfo()
```

    ## R version 4.6.2 (2026-10-29)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 26.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] tidyr_1.3.2      dplyr_1.2.1      purrr_1.2.2      stringr_1.6.0   
    ## [5] rprojroot_2.1.1  knitr_1.52       tibble_3.3.1     BiocStyle_2.41.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] jsonlite_2.0.0      compiler_4.6.2      BiocManager_1.30.27
    ##  [4] tidyselect_1.2.1    jquerylib_0.1.4     systemfonts_1.3.2  
    ##  [7] textshaping_1.0.5   yaml_2.3.12         fastmap_1.2.0      
    ## [10] R6_2.6.1            generics_0.1.4      htmlwidgets_1.6.4  
    ## [13] bookdown_0.48       desc_1.4.3          bslib_0.12.0       
    ## [16] pillar_1.11.1       rlang_1.3.0         cachem_1.1.0       
    ## [19] stringi_1.8.9       xfun_0.62           fs_2.1.0           
    ## [22] sass_0.4.10         otel_0.3.0          cli_3.6.6          
    ## [25] pkgdown_2.2.1       withr_3.0.3         magrittr_2.0.5     
    ## [28] digest_0.6.39       lifecycle_1.0.5     vctrs_0.7.3        
    ## [31] evaluate_1.0.5      glue_1.8.1          ragg_1.5.2         
    ## [34] rmarkdown_2.32      tools_4.6.2         pkgconfig_2.0.3    
    ## [37] htmltools_0.5.9
