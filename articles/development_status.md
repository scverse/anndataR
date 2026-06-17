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
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L76) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L23) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L80) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L146) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L205) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L73) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L208) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L157) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L241) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L117) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L244) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L189) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L101) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L33) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L105) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L211) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L155) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L53) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L159) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L240) |
| raw |  |  |  |  |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L277) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L280) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L270) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L223) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L95) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L226) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L173) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L259) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L123) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L262) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L200) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L128) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L43) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L132) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L226) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L180) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L63) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L184) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L255) |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L51) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L16) | [✅](https://github.com/scverse/anndataR/blob/main/R/HDF5AnnData.R#L55) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-HDF5AnnData.R#L135) |

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
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L57) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L84) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L72) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L254) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L90) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L74) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L93) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L112) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L120) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L94) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L127) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L235) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L164) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L169) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L167) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L200) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L204) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L181) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L207) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L214) |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L246) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L89) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L249) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L132) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L105) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L79) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L108) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L122) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L142) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L96) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L149) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L241) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L184) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L175) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L187) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L207) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L225) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L187) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L228) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L221) |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L36) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L70) | [✅](https://github.com/scverse/anndataR/blob/main/R/ReticulateAnnData.R#L39) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ReticulateAnnData.R#L106) |

### `ZarrAnnData`

| Slot | Getter | Getter test | Setter | Setter test |
|:---|:--:|:--:|:--:|:--:|
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L74) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L23) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L79) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L147) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L203) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L71) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L207) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L160) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L239) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L115) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L243) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L196) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L99) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L104) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L222) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L153) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L51) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L158) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L253) |
| raw |  |  |  |  |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L265) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L269) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L284) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L221) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L93) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L225) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L178) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L252) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L121) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L256) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L209) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L126) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L131) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L237) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L178) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L61) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L183) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L269) |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L49) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L16) | [✅](https://github.com/scverse/anndataR/blob/main/R/ZarrAnnData.R#L54) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-ZarrAnnData.R#L134) |

## Conversion

These tables show the implementation status of conversion between
`AnnData` and other objects.

### `SingleCellExperiment`

| Slot | From | From test | To | To test |
|:---|:--:|:--:|:--:|:--:|
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L314) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L66) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L209) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L64) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L266) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L24) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L239) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L22) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L265) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L18) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L240) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L18) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L327) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L176) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L347) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L202) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L391) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L90) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L249) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L123) |
| raw |  |  |  |  |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L439) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L145) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L257) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L171) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L290) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L45) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L244) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L43) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L289) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L20) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L245) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L16) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L353) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L183) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L348) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L210) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L415) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_SingleCellExperiment.R#L118) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L253) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_SingleCellExperiment.R#L147) |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/from_SingleCellExperiment.R#L142) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/as_SingleCellExperiment.R#L208) |  |

### `Seurat`

| Slot | From | From test | To | To test |
|:---|:--:|:--:|:--:|:--:|
| layers | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L241) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L73) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L229) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L74) |
| obs | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L189) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L35) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L221) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L33) |
| obs_names | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L188) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L31) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L252) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L29) |
| obsm | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L270) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L124) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L511) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L127) |
| obsp | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L340) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L170) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L268) |  |
| raw |  |  |  |  |
| uns | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L413) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L97) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L289) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L99) |
| var | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L216) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L56) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L237) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L55) |
| var_names | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L215) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L29) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L253) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L27) |
| varm | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L300) | [✅](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-from_Seurat.R#L131) | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L512) | [🚧](https://github.com/scverse/anndataR/blob/main/tests/testthat/test-as_Seurat.R#L135) |
| varp | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L376) |  |  |  |
| X | [✅](https://github.com/scverse/anndataR/blob/main/R/from_Seurat.R#L138) |  | [✅](https://github.com/scverse/anndataR/blob/main/R/as_Seurat.R#L228) |  |

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

    ## R version 4.6.0 (2026-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
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
    ## [5] rprojroot_2.1.1  knitr_1.51       tibble_3.3.1     BiocStyle_2.41.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] jsonlite_2.0.0      compiler_4.6.0      BiocManager_1.30.27
    ##  [4] tidyselect_1.2.1    jquerylib_0.1.4     systemfonts_1.3.2  
    ##  [7] textshaping_1.0.5   yaml_2.3.12         fastmap_1.2.0      
    ## [10] R6_2.6.1            generics_0.1.4      htmlwidgets_1.6.4  
    ## [13] bookdown_0.47       desc_1.4.3          bslib_0.11.0       
    ## [16] pillar_1.11.1       rlang_1.2.0         cachem_1.1.0       
    ## [19] stringi_1.8.7       xfun_0.58           fs_2.1.0           
    ## [22] sass_0.4.10         otel_0.2.0          cli_3.6.6          
    ## [25] pkgdown_2.2.0       withr_3.0.2         magrittr_2.0.5     
    ## [28] digest_0.6.39       lifecycle_1.0.5     vctrs_0.7.3        
    ## [31] evaluate_1.0.5      glue_1.8.1          ragg_1.5.2         
    ## [34] rmarkdown_2.31      tools_4.6.0         pkgconfig_2.0.3    
    ## [37] htmltools_0.5.9
