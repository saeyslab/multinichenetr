Multi-Sample Multi-condition Cell-Cell Communication Analysis via
NicheNet: MIS-C application: All-vs-All
================
Robin Browaeys
2022-03-16

<!-- github markdown built using 
rmarkdown::render("vignettes/detailed_analysis_steps_MISC.Rmd", output_format = "github_document")
-->

In this vignette, you can learn how to perform an all-vs-all
MultiNicheNet analysis. In this vignette, we start from one
SingleCellExperiment object containing cells from both sender and
receiver cell types and from different patients.

A MultiNicheNet analysis can be performed if you have multi-sample,
multi-group single-cell data. MultiNicheNet will look for cell-cell
communication between the cell types in your data for each sample, and
compare the cell-cell communication patterns between the groups of
interest. Therefore, the absolute minimum of meta data you need to have,
are following columns indicating for each cell: the **group**,
**sample** and **cell type**.

As example expression data of interacting cells, we will here use
scRNAseq data of immune cells in MIS-C patients and healthy siblings
from this paper of Hoste et al.: [TIM3+ TRBV11-2 T cells and IFNγ
signature in patrolling monocytes and CD16+ NK cells delineate
MIS-C](https://rupress.org/jem/article/219/2/e20211381/212918/TIM3-TRBV11-2-T-cells-and-IFN-signature-in)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6362434.svg)](https://doi.org/10.5281/zenodo.6362434)
. MIS-C (multisystem inflammatory syndrome in children) is a novel rare
immunodysregulation syndrome that can arise after SARS-CoV-2 infection
in children. We will use NicheNet to explore immune cell crosstalk
enriched in MIS-C compared to healthy siblings.

In this vignette, we will prepare the data and analysis parameters, and
then perform the MultiNicheNet analysis.

The different steps of the MultiNicheNet analysis are the following:

-   0.  Preparation of the analysis: load packages, NicheNet LR network
        & ligand-target matrix, single-cell expression data, and define
        main settings of the MultiNicheNet analysis

-   1.  Extract cell type abundance and expression information from
        receiver and sender cell types, and link this expression
        information for ligands of the sender cell types to the
        corresponding receptors of the receiver cell types

-   2.  Perform genome-wide differential expression analysis of receiver
        and sender cell types to define DE genes between the conditions
        of interest. Based on this analysis, we can define the
        logFC/p-value of ligands in senders and receptors in receivers,
        and define the set of affected target genes in the receiver.

-   3.  Predict NicheNet ligand activities and NicheNet ligand-target
        links based on these differential expression results

-   4.  Use the information collected above to prioritize all
        sender-ligand—receiver-receptor pairs.

-   5.  Calculate correlation in expression between ligand-receptor
        pairs and their predicted target genes

In this vignette, we will demonstrate all these steps in detail.

After the MultiNicheNet analysis is done, we will explore the output of
the analysis with different ways of visualization.

# Step 0: Preparation of the analysis: load packages, NicheNet LR network & ligand-target matrix, single-cell expression data

## Step 0.1: Load required packages and NicheNet ligand-receptor network and ligand-target matrix

``` r
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
```

The Nichenet v2 networks and matrices for both mouse and human can be
downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5884439.svg)](https://doi.org/10.5281/zenodo.5884439).

We will read these object in for human because our expression data is of
human patients. Gene names are here made syntactically valid via
`make.names()` to avoid the loss of genes (eg H2-M3) in downstream
visualizations.

``` r
organism = "human"
if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/5884439/files/lr_network_human_21122021.rds"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/5884439/files/ligand_target_matrix_nsga2r_final.rds"))
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/5884439/files/lr_network_mouse_21122021.rds"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/5884439/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
}
```

## Step 0.2: Read in SingleCellExperiment Objects

In this vignette, sender and receiver cell types are in the same
SingleCellExperiment object, which we will load here. In this vignette,
we will load in a subset of the scRNAseq data of the MIS-C
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6362434.svg)](https://doi.org/10.5281/zenodo.6362434).
This subset only contains the cell types that were found to be most
interesting related to MIS-C according to Hoste et al. 

If you start from a Seurat object, you can convert it easily to a
SingleCellExperiment via
`sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")`.

Because the NicheNet 2.0. networks are in the most recent version of the
official gene symbols, we will make sure that the gene symbols used in
the expression data are also updated (= converted from their “aliases”
to official gene symbols). Afterwards, we will make them again
syntactically valid.

``` r
sce = readRDS(url("https://zenodo.org/record/6362434/files/sce_misc_subset.rds"))
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
## [1] "there are provided symbols that are not in the alias annotation table: "
##    [1] "AL627309.1"  "AL669831.5"  "AL645608.8"  "AL645608.2"  "AL162741.1"  "AL391244.3"  "AL645728.1"  "AL691432.2"  "FO704657.1"  "AL109917.1"  "AL590822.2"  "AL139246.5"  "AL139246.3"  "AL513320.1"  "AL805961.2"  "AL365255.1"  "AL031848.2"  "Z98884.1"    "AL034417.4"  "AL034417.3"  "AL096855.1"  "AL139424.2"  "AL354956.1"  "AL139423.1"  "AL109811.3"  "AL109811.1" 
##   [27] "AL109811.2"  "AC243836.1"  "AC254633.1"  "AL031283.1"  "AL121992.3"  "AL121992.1"  "AL450998.2"  "AL137802.2"  "AL021920.3"  "BX284668.2"  "BX284668.5"  "BX284668.6"  "AL035413.1"  "AL031727.2"  "AL020998.1"  "AL031005.1"  "AL031728.1"  "AL031281.2"  "AL445686.2"  "AL445471.1"  "AL445471.2"  "AL031432.5"  "AL031432.3"  "AL606491.1"  "AL031280.1"  "AL020996.1" 
##   [53] "AL033528.2"  "AL391650.1"  "AL512408.1"  "BX293535.1"  "AL512288.1"  "AL353622.1"  "AL353622.2"  "AL360012.1"  "AC114488.2"  "AC114488.3"  "AL136115.2"  "AL445248.1"  "AL049795.1"  "AL662907.3"  "AC114490.3"  "AC114490.2"  "AC004865.2"  "AL591845.1"  "AL139260.1"  "AL033527.3"  "AL033527.5"  "AL050341.2"  "AL603839.3"  "AC098484.3"  "AL512353.1"  "AL139289.2" 
##   [79] "AL139289.1"  "AL451062.3"  "AL357079.1"  "AL139220.2"  "AL592166.1"  "AC104170.1"  "AL050343.1"  "AL445685.1"  "AL606760.2"  "AL606760.3"  "AL606760.1"  "AC119428.2"  "AL161644.1"  "AC105277.1"  "AC099791.2"  "AL357078.1"  "AL583808.1"  "AC118549.1"  "AC103591.3"  "AL359504.2"  "AL078459.1"  "AL049597.2"  "AC099063.4"  "AC099063.1"  "AC099568.2"  "AC104836.1" 
##  [105] "AC105942.1"  "AC092802.1"  "AC118553.1"  "AC118553.2"  "AC104506.1"  "AC093157.2"  "AC093157.1"  "AL109741.1"  "AL390036.1"  "AL359258.2"  "AL356488.3"  "AL356488.2"  "SARS"        "AC000032.1"  "AL450468.2"  "AL160006.1"  "AL355488.1"  "AL365361.1"  "AL360270.1"  "AL360270.3"  "AL355816.2"  "AL354760.1"  "AL603832.1"  "AL357055.3"  "AL137856.1"  "AL390066.1" 
##  [131] "AL445231.1"  "AL157902.1"  "AL359915.2"  "AC244453.2"  "AC244021.1"  "AL592494.3"  "AC239800.2"  "AC239800.3"  "AC245595.1"  "AC246785.3"  "FP700111.1"  "AC245014.3"  "AC245014.1"  "AC243547.3"  "AC243547.2"  "NOTCH2NL"    "AC239799.2"  "AC242426.2"  "AC239809.3"  "AC245297.3"  "AC245297.2"  "AC239868.3"  "AC239868.2"  "AL590133.2"  "AL391069.2"  "AL391069.3" 
##  [157] "AL391335.1"  "AL589765.7"  "AL450992.2"  "AL162258.2"  "AL358472.5"  "AL358472.4"  "AL358472.3"  "AL358472.2"  "AL451085.1"  "AL451085.2"  "AC234582.1"  "AC234582.2"  "AL353807.2"  "AL355388.2"  "AL139412.1"  "AL590560.1"  "AL590560.2"  "AL121987.2"  "AL139011.1"  "AL138930.1"  "AL121985.1"  "AL591806.3"  "AL590714.1"  "AL592295.4"  "AL359541.1"  "AL596325.1" 
##  [183] "AL356441.1"  "AL451074.2"  "AL451050.2"  "AL359962.2"  "AL031733.2"  "AL021068.1"  "AL031599.1"  "AL645568.1"  "AL121983.1"  "AL121983.2"  "Z99127.1"    "AL513329.1"  "AL359265.3"  "AL137796.1"  "AL449106.1"  "AL353708.3"  "AL353708.1"  "AL162431.2"  "AL445228.2"  "AL078644.2"  "AL133383.1"  "AL390957.1"  "AL136987.1"  "AL136454.1"  "AL157402.2"  "AC096677.1" 
##  [209] "AC099676.1"  "AC098934.4"  "AL512306.3"  "AL512306.2"  "AC244034.2"  "AL137789.1"  "AL031316.1"  "AL606468.1"  "AL451060.1"  "AL360091.3"  "AC092803.2"  "AL590648.3"  "AC096642.1"  "AL513314.2"  "AL592148.3"  "AL392172.1"  "AL359979.2"  "AC138393.3"  "AC092809.4"  "AC092811.2"  "AL591895.1"  "AL512343.2"  "AL353593.1"  "AL353593.2"  "AL670729.1"  "AL117350.1" 
##  [235] "AL121990.1"  "AL445524.1"  "AL355472.3"  "AL355472.4"  "AL355472.1"  "AL160408.2"  "AL391832.2"  "AL391832.3"  "AL359921.2"  "AL359921.1"  "BX323046.1"  "AL356512.1"  "AL591848.3"  "AL390728.6"  "AL390728.5"  "AC138089.1"  "AC098483.1"  "AC092159.2"  "AC093390.2"  "AC231981.1"  "AC108488.1"  "AC017076.1"  "AC010904.2"  "AC092580.2"  "AC073195.1"  "AC010969.2" 
##  [261] "AC104794.2"  "AC007240.1"  "AC007249.2"  "AC062028.1"  "AC079148.1"  "AC013400.1"  "AC079145.1"  "AC098828.2"  "AC012065.3"  "AC012065.4"  "AC018742.1"  "AC012073.1"  "AC012074.1"  "AC104699.1"  "AC013472.3"  "AC013403.2"  "AC074117.1"  "AC093690.1"  "AC104695.2"  "AC104695.3"  "AC092164.1"  "AC016907.2"  "AL121652.1"  "AL121658.1"  "AL133245.1"  "AC007378.1" 
##  [287] "AC006369.1"  "AC074366.1"  "AC019171.1"  "AC007388.1"  "AC083949.1"  "AC010883.1"  "AC019129.2"  "AC018682.1"  "AC016722.3"  "AC073283.1"  "AC079807.1"  "AC092650.1"  "AC093635.1"  "AC008280.3"  "AC093110.1"  "AC012358.3"  "AC015982.1"  "AC007743.1"  "AC016727.1"  "AC107081.2"  "AC009501.1"  "AC012368.1"  "AC008074.3"  "AC007365.1"  "AC017083.3"  "AC017083.2" 
##  [313] "AC017083.1"  "AC022201.2"  "AC007881.3"  "AC007878.1"  "AC073046.1"  "AC073263.1"  "AC073263.2"  "AC006030.1"  "AC005041.1"  "AC019069.1"  "AC104135.1"  "AC005034.3"  "AC009309.1"  "AC015971.1"  "AC133644.2"  "AC104134.1"  "AC062029.1"  "AC103563.7"  "AC092835.1"  "AC009237.14" "AC021188.1"  "AC092683.1"  "AC109826.1"  "AC079447.1"  "AC092587.1"  "AC018690.1" 
##  [339] "AC092667.1"  "AC016738.1"  "AC104655.1"  "AC012360.1"  "AC012360.3"  "AC010978.1"  "AC013271.1"  "AC108463.3"  "AC068491.3"  "AC017002.3"  "AC017002.1"  "AC079922.2"  "AC016683.1"  "AC016745.2"  "AC104653.1"  "AC110769.2"  "AC009312.1"  "AC009404.1"  "AC093901.1"  "AC012447.1"  "AC068282.1"  "AC012306.2"  "AC073869.3"  "AC011893.1"  "AC114763.1"  "AC092620.1" 
##  [365] "AC007364.1"  "AC091488.1"  "AC009506.1"  "AC009961.1"  "AC009299.3"  "AC008063.1"  "AC009495.3"  "AC007405.3"  "AC078883.2"  "AC078883.3"  "AC078883.1"  "AC010894.2"  "AC010894.4"  "AC093459.1"  "AC079305.1"  "AC073834.1"  "AC009948.1"  "AC009948.4"  "AC009948.3"  "AC010680.3"  "AC010680.2"  "AC010680.1"  "AC068196.1"  "AC009962.1"  "AC096667.1"  "AC017071.1" 
##  [391] "AC013468.1"  "AC108047.1"  "AC005540.1"  "AC064834.1"  "AC114760.2"  "AC020571.1"  "AC013264.1"  "AC010746.1"  "AC011997.1"  "AC007163.1"  "AC005037.1"  "AC064836.3"  "AC007383.2"  "AC007383.3"  "AC008269.1"  "AC007879.2"  "AC009226.1"  "AC007038.2"  "AC007038.1"  "AC079610.2"  "AC079610.1"  "AC068051.1"  "AC016708.1"  "AC098820.2"  "AC021016.3"  "AC012510.1" 
##  [417] "AC009974.1"  "AC097468.3"  "AC068946.1"  "AC053503.2"  "AC079834.2"  "AC097461.1"  "AC097662.1"  "AC009950.1"  "AC012507.1"  "AC012507.4"  "AC073254.1"  "AC105760.1"  "AC105760.2"  "AC112721.2"  "AC104667.1"  "AC012485.3"  "AC016999.1"  "AC062017.1"  "SEPT2"       "AC114730.1"  "AC131097.4"  "AC093642.1"  "AC024060.1"  "AC018816.1"  "AC026202.2"  "AC026202.3" 
##  [443] "AC026191.1"  "AC018809.2"  "AC034198.2"  "AC093495.1"  "AC090948.1"  "AC090948.2"  "AC090948.3"  "AC116096.1"  "AC112220.4"  "AC123023.1"  "AC097359.3"  "AC097359.2"  "AC092053.3"  "AC092053.2"  "AC099541.1"  "AC006059.1"  "AC099329.3"  "AC104184.1"  "AC104187.1"  "AC099669.1"  "AC124045.1"  "AC098613.1"  "AC099778.1"  "AC134772.1"  "AC137630.2"  "AC137630.1" 
##  [469] "AC137630.3"  "QARS"        "AC105935.1"  "U73166.1"    "AC096920.1"  "AC006252.1"  "AC006254.1"  "AC096887.1"  "AC012467.2"  "AC135507.1"  "AC012557.2"  "AC109587.1"  "AC097634.3"  "AC097634.1"  "AC104435.2"  "AC107204.1"  "AC069222.1"  "AC128688.2"  "AC074044.1"  "AC078785.1"  "AC093010.2"  "AC073352.1"  "AC073352.2"  "AC092910.3"  "AC083798.2"  "AC092903.2" 
##  [495] "AC079848.1"  "AC112484.3"  "AC112484.1"  "AC108673.2"  "AC108673.3"  "AC107027.1"  "AC117382.2"  "AC096992.2"  "AC097103.2"  "AC026304.1"  "AC020636.1"  "AC108718.1"  "AC084036.1"  "AC106707.1"  "AC080013.6"  "AC080013.1"  "AC080013.5"  "AC069224.1"  "AC112491.1"  "AC078802.1"  "AC078795.1"  "AC078795.2"  "AC008040.1"  "AC008040.5"  "AC007823.1"  "AC007620.2" 
##  [521] "AC125618.1"  "AC069431.1"  "AC092953.2"  "AC131235.3"  "AC068631.1"  "AC112907.3"  "AC069213.1"  "AC233280.1"  "AC055764.2"  "AC055764.1"  "AC107464.1"  "AC107464.3"  "AC139887.4"  "AC139887.2"  "AC139887.1"  "AC092535.4"  "AC147067.1"  "AC016773.1"  "AL158068.2"  "AC105415.1"  "AC093323.1"  "AC097382.2"  "AC104825.1"  "AC105345.1"  "AC099550.1"  "AC116651.1" 
##  [547] "AC006160.1"  "AC133961.1"  "AC106047.1"  "AC104078.1"  "AC108471.2"  "AC095057.2"  "AC111006.1"  "AC096734.2"  "AC096586.1"  "AC107068.1"  "AC027271.1"  "AC110792.3"  "AC068620.1"  "AC068620.2"  "AC104806.2"  "AC053527.1"  "AC093677.2"  "AC098818.2"  "AC124016.2"  "AC124016.1"  "AC093827.5"  "AC093827.4"  "AC083829.2"  "AC019077.1"  "AC114811.2"  "AC019131.2" 
##  [573] "AP002026.1"  "AC097460.1"  "AP001816.1"  "AC098487.1"  "AF213884.3"  "AC018797.2"  "AC004069.1"  "AC109361.1"  "AC096564.1"  "AC096564.2"  "AC004067.1"  "AC126283.1"  "AC109347.1"  "AC106864.1"  "AC110079.2"  "AC108866.1"  "AC108062.1"  "AC096711.3"  "AC097376.2"  "AC112236.2"  "AC096733.2"  "AC097504.2"  "AC139720.1"  "AC104596.1"  "AC093835.1"  "AC095055.1" 
##  [599] "AC097375.1"  "AC106882.1"  "AC105285.1"  "AC097534.1"  "AC097534.2"  "AC019163.1"  "AC078881.1"  "AC107214.1"  "AC084871.2"  "AC106897.1"  "AC097652.1"  "AC021087.1"  "AC106772.1"  "AC026740.1"  "AC116351.1"  "AC026412.3"  "AC091965.1"  "AC091965.4"  "AC034229.4"  "AC012640.2"  "AC010491.1"  "AC022113.1"  "AC022113.2"  "AC026785.2"  "AC025181.2"  "AC034231.1" 
##  [625] "AC026801.2"  "AC137810.1"  "AC008957.1"  "AC026741.1"  "AC008875.1"  "AC025171.3"  "AC025171.5"  "AC025171.2"  "AC025171.4"  "AC114947.2"  "AC114956.1"  "AC093297.2"  "AC008966.1"  "AC008914.1"  "AC008892.1"  "AC016596.1"  "AC008937.1"  "AC008937.3"  "AC092343.1"  "AC104113.1"  "AC092354.2"  "AC025442.1"  "AC024581.1"  "AC010280.2"  "AC146944.4"  "AC035140.1" 
##  [651] "AC008972.2"  "AC099522.2"  "AC010501.1"  "AC116337.3"  "AC008897.2"  "AC010245.2"  "AC026725.1"  "AC012636.1"  "AC010260.1"  "AC008771.1"  "AC008434.1"  "AC104118.1"  "AC018754.1"  "AC123595.1"  "AC008840.1"  "AC020900.1"  "AC008906.1"  "AC008906.2"  "AC009126.1"  "AC008522.1"  "AC021086.1"  "AC008467.1"  "AC008494.3"  "AC010226.1"  "AC008549.2"  "AC034236.2" 
##  [677] "AC008629.1"  "AC093535.1"  "AC011416.3"  "AC116366.3"  "AC116366.2"  "AC116366.1"  "AC004775.1"  "AC010240.3"  "AC008608.2"  "AC104109.3"  "AC104109.2"  "AC109454.2"  "AC026691.1"  "AC106791.1"  "AC113382.1"  "AC011405.1"  "AC135457.1"  "AC008667.1"  "AC008438.1"  "AC005618.1"  "AC008781.1"  "AC091959.3"  "AC091948.1"  "AC131025.1"  "AC008641.1"  "AC011337.1" 
##  [703] "AC011374.2"  "AC091982.3"  "AC009185.1"  "AC008676.1"  "AC134043.1"  "AC008691.1"  "AC113414.1"  "AC034199.1"  "AC022217.3"  "AC008429.1"  "AC008429.3"  "AC008378.1"  "AC139491.1"  "AC139491.7"  "AC139795.3"  "AC139795.2"  "AC106795.2"  "AC136604.2"  "AC136604.3"  "AC008393.1"  "AC008610.1"  "AC008443.5"  "AC008443.6"  "AC138035.1"  "AL357054.4"  "AL357054.2" 
##  [729] "AL133351.1"  "AL031963.3"  "AL138831.1"  "AL162718.1"  "AL359643.3"  "AL136307.1"  "AL024498.1"  "AL157373.2"  "AL008729.1"  "AL008729.2"  "AL441883.1"  "AL109914.1"  "AL138720.1"  "AL136162.1"  "AL137003.1"  "AL137003.2"  "AL138724.1"  "AL031775.1"  "U91328.2"    "U91328.1"    "HIST1H2BC"   "AL353759.1"  "HIST1H2BE"   "AL031777.3"  "HIST1H2BF"   "HIST1H2BG"  
##  [755] "HIST1H2BI"   "AL121936.1"  "AL513548.3"  "AL513548.1"  "AL009179.1"  "AL121944.1"  "AL358933.1"  "AL049543.1"  "AL645929.2"  "AL662797.1"  "AL662844.4"  "AL645933.2"  "AL662899.2"  "AL662884.4"  "AL662796.1"  "AL669918.1"  "AL645941.3"  "AL645940.1"  "AL451165.2"  "Z84485.1"    "AL353597.1"  "AL121574.1"  "AL365205.1"  "AL513008.1"  "AL035587.1"  "AL133375.1" 
##  [781] "AL355802.2"  "AL096865.1"  "AL035701.1"  "AL355353.1"  "AL021368.1"  "AL021368.2"  "AL135905.1"  "AL354719.2"  "AL391807.1"  "AL365226.2"  "AL365226.1"  "AC019205.1"  "AL121972.1"  "AL590428.1"  "AL080250.1"  "AL359715.3"  "AL359715.1"  "AL359715.2"  "AL139274.2"  "AL049697.1"  "AL353135.1"  "AL096678.1"  "AL121787.1"  "AL132996.1"  "AL589740.1"  "AL513550.1" 
##  [807] "AL137784.2"  "AL133338.1"  "AL133406.3"  "AL133406.2"  "AL024507.2"  "AL390208.1"  "AL359711.2"  "AL512303.1"  "AC002464.1"  "AL080317.1"  "AL080317.2"  "AL080317.3"  "Z97989.1"    "AL050331.2"  "AL365275.1"  "AL096711.2"  "AL590006.1"  "AL355581.1"  "AL353596.1"  "AL024508.2"  "AL356234.2"  "AL357060.1"  "AL590617.2"  "AL031772.1"  "AL158850.1"  "AL035446.1" 
##  [833] "AL023806.1"  "AL356599.1"  "AL031056.1"  "AL358852.1"  "AL355312.4"  "AL080276.2"  "AL590867.1"  "AL596202.1"  "AL355297.3"  "AL355297.4"  "AL391863.1"  "AL035634.1"  "AL627422.2"  "AL035530.2"  "AL356417.3"  "AL356417.1"  "AL139393.2"  "AL022069.1"  "AL159163.1"  "Z94721.1"    "AL121935.1"  "AL121935.2"  "BX322234.1"  "AL109910.2"  "AC093627.5"  "AC093627.4" 
##  [859] "AC147651.1"  "AC147651.4"  "AC073957.3"  "AC073957.2"  "AC091729.3"  "AC102953.2"  "AC104129.1"  "AC004906.1"  "AC024028.1"  "AC092171.4"  "AC092171.3"  "AC073343.2"  "AC004982.2"  "AC004948.1"  "AC006042.4"  "AC006042.2"  "AC007029.1"  "AC073333.1"  "AC073332.1"  "AC004130.1"  "AC073072.1"  "AC005082.1"  "MPP6"        "AC005165.1"  "AC004520.1"  "AC004540.1" 
##  [885] "AC004593.3"  "AC007285.1"  "AC007036.1"  "AC005154.6"  "GARS"        "AC018647.2"  "AC007349.3"  "AC006033.2"  "AC011290.1"  "AC072061.1"  "AC004951.1"  "AC004854.2"  "AC073115.2"  "AC020743.1"  "AC074351.1"  "AC118758.3"  "AC092634.3"  "AC073349.1"  "AC068533.4"  "AC008267.5"  "AC027644.3"  "AC073335.2"  "AC006480.2"  "AC211476.2"  "AC004921.1"  "AC005076.1" 
##  [911] "AC003991.2"  "AC002383.1"  "AC007566.1"  "AC002454.1"  "AC002451.1"  "AC004834.1"  "AC004922.1"  "AC073842.2"  "AC073842.1"  "AC092849.1"  "AC069281.1"  "AC105446.1"  "AC006329.1"  "AC093668.3"  "AC005064.1"  "AC007384.1"  "AC005070.3"  "AC073073.2"  "AC007032.1"  "AC004492.1"  "AC004839.1"  "AC002467.1"  "AC005046.1"  "AC006333.2"  "AC073934.1"  "AC090114.2" 
##  [937] "AC018638.7"  "AC078846.1"  "AC073320.1"  "AC087071.1"  "AC007938.2"  "AC007938.3"  "AC016831.7"  "AC016831.1"  "AC016831.5"  "AC058791.1"  "AC008264.2"  "AC009275.1"  "AC083862.2"  "AC009542.1"  "AC024084.1"  "AC078845.1"  "AC083880.1"  "AC004918.1"  "AC093673.1"  "AC004889.1"  "AC074386.1"  "AC005229.4"  "AC073314.1"  "AC004877.1"  "AC092681.3"  "AC005586.1" 
##  [963] "AC073111.5"  "AC021097.1"  "AC006017.1"  "AC093726.2"  "AC144652.1"  "AC009403.1"  "AC005534.1"  "AC005481.1"  "AC011899.2"  "AC011899.3"  "AL732314.6"  "AL732314.4"  "AL672277.1"  "AL683807.1"  "BX890604.1"  "AC073529.1"  "AC131011.1"  "CXorf21"     "AC108879.1"  "AC234772.3"  "AC234772.2"  "AF196972.1"  "AC115618.1"  "AC231533.1"  "AC232271.1"  "AC233976.1" 
##  [989] "AL445472.1"  "AL354793.1"  "AL022157.1"  "AL034397.2"  "AL034397.3"  "AL590764.1"  "AL353804.1"  "Z68871.1"    "AL121601.1"  "AL683813.1"  "AL078639.1"  "AC244197.3" 
##  [ reached getOption("max.print") -- omitted 2569 entries ]
## [1] "they are added to the alias annotation table, so they don't get lost"
## [1] "following are the official gene symbols of input aliases: "
##              symbol            alias
## 1             AARS1             AARS
## 2           ABITRAM          FAM206A
## 3              ACP3             ACPP
## 4          ACTN1-DT        ACTN1-AS1
## 5             ADPRS          ADPRHL2
## 6             ADSS1           ADSSL1
## 7             ADSS2             ADSS
## 8        ANKRD20A2P        ANKRD20A2
## 9        ANKRD20A3P        ANKRD20A3
## 10       ANKRD20A4P        ANKRD20A4
## 11        ANKRD40CL        LINC00483
## 12           ANTKMT          FAM173A
## 13            AOPEP           C9orf3
## 14           ARLNC1        LINC02170
## 15           ATP5MJ          ATP5MPL
## 16           ATP5MK           ATP5MD
## 17         ATPSCKMT          FAM173B
## 18          B3GALT9          B3GNT10
## 19             BBLN          C9orf16
## 20           BMERB1         C16orf45
## 21            BPNT2           IMPAD1
## 22            BRME1         C19orf57
## 23      CARNMT1-AS1      C9orf41-AS1
## 24            CARS1             CARS
## 25        CARS1-AS1         CARS-AS1
## 26             CBY2            SPERT
## 27          CCDC200        LINC00854
## 28      CCDC28A-AS1            GVQW2
## 29             CCN2             CTGF
## 30             CCN3              NOV
## 31             CCNP            CNTD2
## 32            CDIN1         C15orf41
## 33          CENATAC           CCDC84
## 34            CEP20            FOPNL
## 35            CEP43          FGFR1OP
## 36            CERT1         COL4A3BP
## 37         CFAP20DC          C3orf67
## 38          CFAP251            WDR66
## 39          CFAP410          C21orf2
## 40           CFAP91           MAATS1
## 41           CFAP92         KIAA1257
## 42           CIAO2A           FAM96A
## 43           CIAO2B           FAM96B
## 44            CIAO3            NARFL
## 45           CIBAR1           FAM92A
## 46           CIBAR2           FAM92B
## 47            CILK1              ICK
## 48        CLLU1-AS1          CLLU1OS
## 49             COA8           APOPT1
## 50            CRACD         KIAA1211
## 51           CRACDL        KIAA1211L
## 52            CRPPA             ISPD
## 53            CYRIA           FAM49A
## 54            CYRIB           FAM49B
## 55             CZIB         C1orf123
## 56            DARS1             DARS
## 57        DARS1-AS1         DARS-AS1
## 58          DENND10           FAM45A
## 59          DENND11         KIAA1147
## 60          DENND2B              ST5
## 61           DIPK1A           FAM69A
## 62           DIPK1B           FAM69B
## 63           DIPK2A          C3orf58
## 64           DMAC2L            ATP5S
## 65          DNAAF10            WDR92
## 66          DNAAF11            LRRC6
## 67           DNAAF8         C16orf71
## 68           DNAAF9        C20orf194
## 69            DNAI3            WDR63
## 70            DNAI4            WDR78
## 71            DNAI7            CASC1
## 72        DOCK8-AS1          C9orf66
## 73            DOP1A           DOPEY1
## 74            DOP1B           DOPEY2
## 75           DUSP29            DUPD1
## 76          DYNC2I1            WDR60
## 77          DYNC2I2            WDR34
## 78           DYNLT2            TCTE3
## 79          DYNLT2B         TCTEX1D2
## 80           DYNLT4         TCTEX1D4
## 81           DYNLT5         TCTEX1D1
## 82            ECRG4          C2orf40
## 83          ELAPOR1         KIAA1324
## 84          ELAPOR2        KIAA1324L
## 85            EOLA1         CXorf40A
## 86            EOLA2         CXorf40B
## 87            EPRS1             EPRS
## 88            EZHIP          CXorf67
## 89         FAM153CP          FAM153C
## 90          FAM174C         C19orf24
## 91         FAM86C1P          FAM86C1
## 92             FCSK              FUK
## 93           FHIP1A         FAM160A1
## 94        FHIP1A-DT      FAM160A1-DT
## 95           FHIP1B         FAM160A2
## 96           FHIP2A         FAM160B1
## 97           FHIP2B         FAM160B2
## 98           FLACC1         ALS2CR12
## 99      FMC1-LUC7L2   C7orf55-LUC7L2
## 100          GARRE1         KIAA0355
## 101          GASK1A          FAM198A
## 102          GASK1B          FAM198B
## 103      GASK1B-AS1      FAM198B-AS1
## 104           GCNT4        LINC01336
## 105        GDF5-AS1           GDF5OS
## 106            GET1              WRB
## 107            GET3            ASNA1
## 108            GFUS            TSTA3
## 109           GOLM2            CASC4
## 110            H1-0             H1F0
## 111            H1-1         HIST1H1A
## 112           H1-10             H1FX
## 113       H1-10-AS1         H1FX-AS1
## 114            H1-2         HIST1H1C
## 115            H1-3         HIST1H1D
## 116            H1-4         HIST1H1E
## 117            H1-5         HIST1H1B
## 118            H1-6         HIST1H1T
## 119           H2AB1           H2AFB1
## 120          H2AC11        HIST1H2AG
## 121          H2AC12        HIST1H2AH
## 122          H2AC13        HIST1H2AI
## 123          H2AC14        HIST1H2AJ
## 124          H2AC15        HIST1H2AK
## 125          H2AC16        HIST1H2AL
## 126          H2AC17        HIST1H2AM
## 127          H2AC18       HIST2H2AA3
## 128          H2AC19       HIST2H2AA4
## 129          H2AC20        HIST2H2AC
## 130          H2AC21        HIST2H2AB
## 131           H2AC4        HIST1H2AB
## 132           H2AC6        HIST1H2AC
## 133           H2AC8        HIST1H2AE
## 134            H2AJ            H2AFJ
## 135            H2AW         HIST3H2A
## 136            H2AX            H2AFX
## 137           H2AZ1            H2AFZ
## 138           H2AZ2            H2AFV
## 139          H2BC11        HIST1H2BJ
## 140          H2BC12        HIST1H2BK
## 141          H2BC13        HIST1H2BL
## 142          H2BC14        HIST1H2BM
## 143          H2BC15        HIST1H2BN
## 144          H2BC17        HIST1H2BO
## 145          H2BC18        HIST2H2BF
## 146          H2BC21        HIST2H2BE
## 147           H2BC3        HIST1H2BB
## 148           H2BC5        HIST1H2BD
## 149           H2BC9        HIST1H2BH
## 150           H2BU1        HIST3H2BB
## 151            H3-2       HIST2H3PS2
## 152           H3-3A            H3F3A
## 153           H3-3B            H3F3B
## 154            H3-5            H3F3C
## 155            H3C1         HIST1H3A
## 156           H3C10         HIST1H3H
## 157           H3C11         HIST1H3I
## 158           H3C12         HIST1H3J
## 159           H3C13         HIST2H3D
## 160            H3C2         HIST1H3B
## 161            H3C3         HIST1H3C
## 162            H3C6         HIST1H3E
## 163            H3C7         HIST1H3F
## 164            H3C8         HIST1H3G
## 165           H4-16          HIST4H4
## 166            H4C1         HIST1H4A
## 167           H4C11         HIST1H4J
## 168           H4C13         HIST1H4L
## 169           H4C14         HIST2H4A
## 170           H4C15         HIST2H4B
## 171            H4C2         HIST1H4B
## 172            H4C3         HIST1H4C
## 173            H4C4         HIST1H4D
## 174            H4C5         HIST1H4E
## 175            H4C6         HIST1H4F
## 176            H4C8         HIST1H4H
## 177            H4C9         HIST1H4I
## 178           HARS1             HARS
## 179            HEXD            HEXDC
## 180            HROB         C17orf53
## 181       HSDL2-AS1         C9orf147
## 182           IARS1             IARS
## 183           IFTAP         C11orf74
## 184           ILRUN         C6orf106
## 185           IRAG1            MRVI1
## 186       IRAG1-AS1        MRVI1-AS1
## 187           IRAG2             LRMP
## 188         ITPRID2            SSFA2
## 189           KARS1             KARS
## 190           KASH5          CCDC155
## 191          KATNIP         KIAA0556
## 192           KICS2         C12orf66
## 193           KIFBP           KIF1BP
## 194       KRT10-AS1           TMEM99
## 195           LARS1             LARS
## 196           LDAF1          TMEM159
## 197       LINC02481       LINC002481
## 198       LINC02693         C17orf51
## 199       LINC02694         C15orf53
## 200       LINC02870         C10orf91
## 201       LINC02875         C17orf82
## 202       LINC02899          C5orf17
## 203       LINC02901          C6orf99
## 204       LINC02908         C9orf139
## 205       LINC02909         C12orf77
## 206       LINC02910        C20orf197
## 207       LINC02913         C9orf106
## 208       LNCTAM34A        LINC01759
## 209          LRATD2           FAM84B
## 210            LTO1           ORAOV1
## 211           MACIR          C5orf30
## 212       MACROH2A1            H2AFY
## 213       MACROH2A2           H2AFY2
## 214           MAILR        AZIN1-AS1
## 215         MARCHF1           MARCH1
## 216        MARCHF10          MARCH10
## 217         MARCHF2           MARCH2
## 218         MARCHF3           MARCH3
## 219         MARCHF4           MARCH4
## 220         MARCHF5           MARCH5
## 221         MARCHF6           MARCH6
## 222         MARCHF7           MARCH7
## 223         MARCHF8           MARCH8
## 224         MARCHF9           MARCH9
## 225           MEAK7            TLDC1
## 226        METTL25B           RRNAD1
## 227         MICOS10           MINOS1
## 228         MICOS13         C19orf70
## 229          MIDEAS          ELMSAN1
## 230          MINAR1         KIAA1024
## 231       MIR3667HG         C22orf34
## 232        MIR9-1HG          C1orf61
## 233           MIX23           CCDC58
## 234            MMUT              MUT
## 235          MROCKI        LINC01268
## 236           MRTFA             MKL1
## 237           MRTFB             MKL2
## 238          MTARC1            MARC1
## 239          MTARC2            MARC2
## 240            MTLN           SMIM37
## 241          MTRES1         C6orf203
## 242           MTRFR         C12orf65
## 243           MTSS2           MTSS1L
## 244            MYG1         C12orf10
## 245           NARS1             NARS
## 246        NCBP2AS2        NCBP2-AS2
## 247       NDUFV1-DT         C11orf72
## 248          NIBAN1          FAM129A
## 249          NIBAN2          FAM129B
## 250          NIBAN3          FAM129C
## 251        NOPCHAP1         C12orf45
## 252           NTAQ1           WDYHV1
## 253           NUP42            NUPL2
## 254            OBI1           RNF219
## 255        OBI1-AS1       RNF219-AS1
## 256           ODAD1          CCDC114
## 257           ODAD2            ARMC4
## 258           ODAD3          CCDC151
## 259           ODAD4            TTC25
## 260          PABIR1          FAM122A
## 261          PABIR2          FAM122B
## 262          PABIR3          FAM122C
## 263           PACC1          TMEM206
## 264      PALM2AKAP2            AKAP2
## 265      PALM2AKAP2      PALM2-AKAP2
## 266           PALS1             MPP5
## 267       PAXIP1-DT       PAXIP1-AS1
## 268           PEDS1          TMEM189
## 269    PEDS1-UBE2V1   TMEM189-UBE2V1
## 270         PELATON           SMIM25
## 271           PGAP4          TMEM246
## 272           PGAP6           TMEM8A
## 273           PHAF1         C16orf70
## 274      PIK3IP1-DT      PIK3IP1-AS1
## 275       PITX1-AS1          C5orf66
## 276         PLA2G2C       UBXN10-AS1
## 277          PLAAT1           HRASLS
## 278          PLAAT2          HRASLS2
## 279          PLAAT3          PLA2G16
## 280          PLAAT4          RARRES3
## 281          PLAAT5          HRASLS5
## 282         PLEKHG7         C12orf74
## 283         POGLUT2           KDELC1
## 284         POGLUT3           KDELC2
## 285          POLR1F          TWISTNB
## 286          POLR1G           CD3EAP
## 287          POLR1H            ZNRD1
## 288     PPP1R13B-DT        LINC00637
## 289          PRANCR        LINC01481
## 290         PRECSIT        LINC00346
## 291           PRORP         KIAA0391
## 292         PRSS42P           PRSS42
## 293          PRXL2A          FAM213A
## 294          PRXL2B          FAM213B
## 295          PRXL2C            AAED1
## 296        PSME3IP1          FAM192A
## 297        RAB30-DT        RAB30-AS1
## 298            RADX          CXorf57
## 299           RAMAC           RAMMET
## 300           RARS1             RARS
## 301            RBIS          C8orf59
## 302           RELCH         KIAA1468
## 303           RESF1         KIAA1551
## 304            RO60           TROVE2
## 305        RPL34-DT        RPL34-AS1
## 306            RRM2          C2orf48
## 307            RSKR           SGK494
## 308           RUSF1         C16orf58
## 309           SANBR         KIAA1841
## 310           SCAT1        LINC02081
## 311         SEPTIN1            SEPT1
## 312        SEPTIN10           SEPT10
## 313        SEPTIN11           SEPT11
## 314         SEPTIN3            SEPT3
## 315         SEPTIN4            SEPT4
## 316     SEPTIN4-AS1        SEPT4-AS1
## 317         SEPTIN5            SEPT5
## 318         SEPTIN6            SEPT6
## 319         SEPTIN7            SEPT7
## 320      SEPTIN7-DT        SEPT7-AS1
## 321         SEPTIN8            SEPT8
## 322         SEPTIN9            SEPT9
## 323            SHFL         C19orf66
## 324           SHOC1          C9orf84
## 325         SLC49A4            DIRC2
## 326         SLC66A1            PQLC2
## 327         SLC66A2            PQLC1
## 328         SLC66A3            PQLC3
## 329         SMC5-DT         SMC5-AS1
## 330          SMIM43          TMEM155
## 331          SNHG30        LINC02001
## 332          SNHG32          C6orf48
## 333         SPRING1         C12orf49
## 334           SSPOP             SSPO
## 335         STAM-DT         STAM-AS1
## 336          STEEP1          CXorf56
## 337  STIMATE-MUSTN1   TMEM110-MUSTN1
## 338          STING1          TMEM173
## 339        STX17-DT        STX17-AS1
## 340           TAFA1          FAM19A1
## 341           TAFA2          FAM19A2
## 342         TAMALIN            GRASP
## 343           TARS1             TARS
## 344           TARS3           TARSL2
## 345           TASOR          FAM208A
## 346          TASOR2          FAM208B
## 347        TBC1D29P          TBC1D29
## 348   TIMM23B-AGAP6        LINC00843
## 349          TLCD3A           FAM57A
## 350          TLCD3B           FAM57B
## 351           TLCD4           TMEM56
## 352           TLCD5          TMEM136
## 353            TLE5              AES
## 354 TM4SF19-DYNLT2B TM4SF19-TCTEX1D2
## 355        TMCC1-DT        TMCC1-AS1
## 356       TOLLIP-DT       TOLLIP-AS1
## 357        TRAPPC14          C7orf43
## 358          TUBB8B         TUBB8P12
## 359       USP27X-DT       USP27X-AS1
## 360        USP46-DT        USP46-AS1
## 361           VARS1             VARS
## 362           WARS1             WARS
## 363            YAE1           YAE1D1
## 364           YARS1             YARS
## 365           YJU2B          CCDC130
## 366       YTHDF3-DT       YTHDF3-AS1
## 367            ZFTA         C11orf95
## 368       ZNF22-AS1         C10orf25
## 369      ZNF407-AS1        LINC00909
## 370       ZNF516-DT         C18orf65
## 371       ZNF582-DT       ZNF582-AS1
## 372          ZNF875             HKR1
## 373           ZNRD2           SSSCA1
```

In this case study, we want to study differences in cell-cell
communication patterns between MIS-C patients (M), their healthy
siblings (S) and adult patients with severe covid (A). The meta data
columns that indicate this disease status is `MIS.C.AgeTier`.

Cell type annotations are indicated in the `Annotation_v2.0` column, and
the sample is indicated by the `ShortID` column. If your cells are
annotated in multiple hierarchical levels, we recommend using a high
level in the hierarchy. This for 2 reasons: 1) MultiNicheNet focuses on
differential expression and not differential abundance, and 2) there
should be sufficient cells per sample-celltype combination.

We will now also check the number of cells per cell type condition
combination, and the number of patients per condition.

``` r
table(SummarizedExperiment::colData(sce)$Annotation_v2.0, SummarizedExperiment::colData(sce)$ShortID) # cell types vs samples
##                          
##                             A1   A2   A3   A4   M1   M2   M3   M4   M5   M6   M7   S1   S2   S3   S4   S5
##   L_NK_CD56._CD16.        1302  712  462 1273  158   53  273  573   96   47  367  620  355  352 1082  329
##   L_NK_CD56hi               33   59   49   35   14    5   29   56    8    7   52   39  126   50   91  169
##   L_T_Proliferating        104  192   97   39  221   27  155   84   69  154 1525   34   33   53   29   21
##   L_T_TIM3._CD38._HLADR.    27  137  291   36  873  102  426  130   56  320 1193   76   38   50   21   23
##   M_Monocyte_CD14_resting   61  436  771   63 1003  750 1120  379   54   45  190  574  956 1396  493 1398
##   M_Monocyte_CD16            6   14   63    3   21   74    5  224   96    1  196  116  104  329   71  167
```

As you can see, some Celltype-Sample combinations have 0 cells. It is
possible that during DE analysis, some cell types will be removed from
the analysis if there is not enough information to do a DE analysis.
(More info later)

``` r
table(SummarizedExperiment::colData(sce)$Annotation_v2.0, SummarizedExperiment::colData(sce)$MIS.C.AgeTier) # cell types vs conditions
##                          
##                              A    M    S
##   L_NK_CD56._CD16.        3749 1567 2738
##   L_NK_CD56hi              176  171  475
##   L_T_Proliferating        432 2235  170
##   L_T_TIM3._CD38._HLADR.   491 3100  208
##   M_Monocyte_CD14_resting 1331 3541 4817
##   M_Monocyte_CD16           86  617  787
```

## Step 0.3: Prepare settings of the MultiNicheNet cell-cell communication analysis

### Define in which metadata columns we can find the **group**, **sample** and **cell type** IDs

If you would have batch effects or covariates you can correct for, you
can define this here as well.

Important: for categorical covariates and batches, there should be at
least one sample for every group-batch combination. If one of your
groups/conditions lacks a certain level of your batch, you won’t be able
to correct for the batch effect because the model is then not able to
distinguish batch from group/condition effects.

Important: the column names of group, sample, cell type, batches and
covariates should be syntactically valid (`make.names`)

Important: All group, sample, cell type, batch and covariate names
should be syntactically valid as well (`make.names`) (eg through
`SummarizedExperiment::colData(sce)$ShortID = SummarizedExperiment::colData(sce)$ShortID %>% make.names()`)

``` r
sample_id = "ShortID"
group_id = "MIS.C.AgeTier"
celltype_id = "Annotation_v2.0"
covariates = NA
batches = NA
```

Sender and receiver cell types also need to be defined. Both are here
all cell types in the dataset because we are interested in an All-vs-All
analysis.

``` r
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
```

If the user wants it, it is possible to use only a subset of senders and
receivers. Senders and receivers can be entirely different, but also
overlapping, or the same. If you don’t use all the cell types in your
data, we recommend to continue with a subset of your data.

``` r
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]
```

Now we will go to the first real step of the MultiNicheNet analysis

# Step 1: Extract cell type abundance and expression information from receiver and sender cell types, and link this expression information for ligands of the sender cell types to the corresponding receptors of the receiver cell types

Since MultiNicheNet will infer group differences at the sample level for
each cell type (currently via Muscat - pseudobulking + EdgeR), we need
to have sufficient cells per sample of a cell type, and this for both
groups. In the following analysis we will set this minimum number of
cells per cell type per sample at 10 (absolute minimum).

``` r
min_cells = 10
```

Now we will calculate abundance and expression information for each cell
type / sample / group combination with the following functions. In the
output of this function, you can also find some ‘Cell type abundance
diagnostic plots’ that will the users which celltype-sample combinations
will be left out later on for DE calculation because the nr of cells is
lower than de defined minimum defined here above. If too many
celltype-sample combinations don’t pass this threshold, we recommend to
define your cell types in a more general way (use one level higher of
the cell type ontology hierarchy) (eg TH17 CD4T cells –&gt; CD4T cells).

``` r
abundance_expression_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network, batches = batches)
```

First, we will check the cell type abundance diagnostic plots.

### Interpretation of cell type abundance information

The first plot visualizes the number of cells per celltype-sample
combination, and indicates which combinations are removed during the DE
analysis because there are less than `min_cells` in the celltype-sample
combination.

``` r
abundance_expression_info$abund_plot_sample
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
The red dotted line indicates the required minimum of cells as defined
above in `min_cells`. We can see here that some sample-celltype
combinations are left out. For the DE analysis in the next step, only
cell types will be considered if there are at least two samples per
group with a sufficient number of cells.

In a next plot, we will look at differential abundance between the
conditions. This because the pseudobulking approach behind Muscat could
potentially suffer from some biases if there would be huge differences
in abundances of a cell type between different groups. Downstream
results of these cell types should then be considered with some caution.

``` r
abundance_expression_info$abund_plot_group
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
Differential abundance looks quite OK.

If you want to look at the cell numbers behind these plots, you can do
so via the following piece of code

``` r
abundance_expression_info$abundance_data_receiver %>% head()
## # A tibble: 6 x 5
## # Groups:   sample, receiver [6]
##   sample receiver                n_cells_receiver group keep_receiver
##   <chr>  <chr>                              <int> <chr>         <dbl>
## 1 A1     L_NK_CD56._CD16.                    1302 A                 1
## 2 A1     L_NK_CD56hi                           33 A                 1
## 3 A1     L_T_Proliferating                    104 A                 1
## 4 A1     L_T_TIM3._CD38._HLADR.                27 A                 1
## 5 A1     M_Monocyte_CD14_resting               61 A                 1
## 6 A1     M_Monocyte_CD16                        6 A                 0
abundance_expression_info$abundance_data_sender %>% head() # in the case of an all-vs-all analysis: both are the same
## # A tibble: 6 x 5
## # Groups:   sample, sender [6]
##   sample sender                  n_cells_sender group keep_sender
##   <chr>  <chr>                            <int> <chr>       <dbl>
## 1 A1     L_NK_CD56._CD16.                  1302 A               1
## 2 A1     L_NK_CD56hi                         33 A               1
## 3 A1     L_T_Proliferating                  104 A               1
## 4 A1     L_T_TIM3._CD38._HLADR.              27 A               1
## 5 A1     M_Monocyte_CD14_resting             61 A               1
## 6 A1     M_Monocyte_CD16                      6 A               0
```

**Important**: Based on the cell type abundance diagnostics, we
recommend users to change their analysis settings if required (eg
changing cell type annotation level, batches, …), before proceeding with
the rest of the analysis.

### Interpretation of expression information

Previously, we also calculated expression information. With the
following piece of code, you can check the average expression for each
gene per sample (normalized expression value and fraction of expressing
cells with non-zero counts, and logCPM-pseudocounts).

``` r
abundance_expression_info$celltype_info$avg_df %>% head()
## # A tibble: 6 x 4
##   gene       sample average_sample celltype              
##   <chr>      <chr>           <dbl> <fct>                 
## 1 AL627309.1 M1            0       L_T_TIM3._CD38._HLADR.
## 2 AL669831.5 M1            0.0231  L_T_TIM3._CD38._HLADR.
## 3 LINC00115  M1            0.0161  L_T_TIM3._CD38._HLADR.
## 4 FAM41C     M1            0.00367 L_T_TIM3._CD38._HLADR.
## 5 NOC2L      M1            0.212   L_T_TIM3._CD38._HLADR.
## 6 KLHL17     M1            0.00913 L_T_TIM3._CD38._HLADR.
abundance_expression_info$celltype_info$frq_df %>% head()
## # A tibble: 6 x 4
##   gene       sample fraction_sample celltype              
##   <chr>      <chr>            <dbl> <chr>                 
## 1 AL627309.1 M1             0       L_T_TIM3._CD38._HLADR.
## 2 AL669831.5 M1             0.0218  L_T_TIM3._CD38._HLADR.
## 3 LINC00115  M1             0.0149  L_T_TIM3._CD38._HLADR.
## 4 FAM41C     M1             0.00344 L_T_TIM3._CD38._HLADR.
## 5 NOC2L      M1             0.190   L_T_TIM3._CD38._HLADR.
## 6 KLHL17     M1             0.00916 L_T_TIM3._CD38._HLADR.
abundance_expression_info$celltype_info$pb_df %>% head()
## # A tibble: 6 x 4
##   gene       sample pb_sample celltype              
##   <chr>      <chr>      <dbl> <fct>                 
## 1 AL627309.1 A1          0    L_T_TIM3._CD38._HLADR.
## 2 AL669831.5 A1          3.11 L_T_TIM3._CD38._HLADR.
## 3 LINC00115  A1          0    L_T_TIM3._CD38._HLADR.
## 4 FAM41C     A1          0    L_T_TIM3._CD38._HLADR.
## 5 NOC2L      A1          4.58 L_T_TIM3._CD38._HLADR.
## 6 KLHL17     A1          0    L_T_TIM3._CD38._HLADR.
```

Now for the average per group:

``` r
abundance_expression_info$celltype_info$avg_df_group %>% head()
## # A tibble: 6 x 4
## # Groups:   group, celltype [1]
##   group celltype         gene     average_group
##   <chr> <fct>            <chr>            <dbl>
## 1 A     L_NK_CD56._CD16. A1BG           0.0537 
## 2 A     L_NK_CD56._CD16. A1BG.AS1       0.0156 
## 3 A     L_NK_CD56._CD16. A2M            0.00100
## 4 A     L_NK_CD56._CD16. A2M.AS1        0.00794
## 5 A     L_NK_CD56._CD16. A4GALT         0      
## 6 A     L_NK_CD56._CD16. AAAS           0.108
abundance_expression_info$celltype_info$frq_df_group %>% head()
## # A tibble: 6 x 4
## # Groups:   group, celltype [1]
##   group celltype         gene     fraction_group
##   <chr> <chr>            <chr>             <dbl>
## 1 A     L_NK_CD56._CD16. A1BG           0.0356  
## 2 A     L_NK_CD56._CD16. A1BG.AS1       0.0103  
## 3 A     L_NK_CD56._CD16. A2M            0.000580
## 4 A     L_NK_CD56._CD16. A2M.AS1        0.00533 
## 5 A     L_NK_CD56._CD16. A4GALT         0       
## 6 A     L_NK_CD56._CD16. AAAS           0.0662
abundance_expression_info$celltype_info$pb_df_group %>% head()
## # A tibble: 6 x 4
## # Groups:   group, celltype [1]
##   group celltype         gene     pb_group
##   <chr> <fct>            <chr>       <dbl>
## 1 A     L_NK_CD56._CD16. A1BG        3.93 
## 2 A     L_NK_CD56._CD16. A1BG.AS1    2.32 
## 3 A     L_NK_CD56._CD16. A2M         0.267
## 4 A     L_NK_CD56._CD16. A2M.AS1     1.39 
## 5 A     L_NK_CD56._CD16. A4GALT      0    
## 6 A     L_NK_CD56._CD16. AAAS        4.83
```

In the last part of this step, we combined this information for each
ligand-receptor pair combination for each sender-receiver combination.
The output of this can be seen as well:

For sample-based:

``` r
abundance_expression_info$sender_receiver_info$avg_df %>% head()
## # A tibble: 6 x 8
##   sample sender                  receiver         ligand receptor avg_ligand avg_receptor ligand_receptor_prod
##   <chr>  <fct>                   <fct>            <chr>  <chr>         <dbl>        <dbl>                <dbl>
## 1 M4     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A8 CD68           6.48         2.81                 18.2
## 2 A1     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A8 CD68           6.33         2.82                 17.9
## 3 M4     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A9 CD68           6.00         2.81                 16.9
## 4 M4     L_NK_CD56._CD16.        M_Monocyte_CD16  B2M    LILRB1         5.07         3.01                 15.3
## 5 A1     M_Monocyte_CD14_resting L_NK_CD56._CD16. S100A8 ITGB2          6.33         2.41                 15.2
## 6 M5     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A8 CD68           5.85         2.56                 15.0
abundance_expression_info$sender_receiver_info$frq_df %>% head()
## # A tibble: 6 x 8
##   sample sender                 receiver        ligand receptor fraction_ligand fraction_receptor ligand_receptor_fraction_prod
##   <chr>  <chr>                  <chr>           <chr>  <chr>              <dbl>             <dbl>                         <dbl>
## 1 M1     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 HLA.C  LILRB1                 1                 1                             1
## 2 M1     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 HLA.C  LILRB2                 1                 1                             1
## 3 M1     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 B2M    LILRB1                 1                 1                             1
## 4 M1     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 RPS19  C5AR1                  1                 1                             1
## 5 M2     L_T_TIM3._CD38._HLADR. L_NK_CD56hi     HLA.A  KLRD1                  1                 1                             1
## 6 M2     L_T_TIM3._CD38._HLADR. L_NK_CD56hi     HLA.B  KLRD1                  1                 1                             1
abundance_expression_info$sender_receiver_info$pb_df %>% head()
## # A tibble: 6 x 8
##   sample sender                  receiver                ligand receptor pb_ligand pb_receptor ligand_receptor_pb_prod
##   <chr>  <fct>                   <fct>                   <chr>  <chr>        <dbl>       <dbl>                   <dbl>
## 1 M4     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A8 CD68          16.1        11.1                    179.
## 2 M4     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A9 CD68          15.8        11.1                    174.
## 3 M4     M_Monocyte_CD14_resting L_NK_CD56._CD16.        S100A8 ITGB2         16.1        10.1                    164.
## 4 M4     M_Monocyte_CD14_resting M_Monocyte_CD14_resting S100A8 ITGB2         16.1        10.1                    163.
## 5 A1     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A8 CD68          15.8        10.3                    163.
## 6 M4     M_Monocyte_CD16         M_Monocyte_CD16         B2M    LILRB1        14.2        11.5                    162.
```

For group-based:

``` r
abundance_expression_info$sender_receiver_info$avg_df_group %>% head()
## # A tibble: 6 x 8
## # Groups:   group, sender [2]
##   group sender                  receiver         ligand receptor avg_ligand_group avg_receptor_group ligand_receptor_prod_group
##   <chr> <fct>                   <fct>            <chr>  <chr>               <dbl>              <dbl>                      <dbl>
## 1 M     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A8 CD68                 5.73               2.35                       13.5
## 2 M     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A9 CD68                 5.61               2.35                       13.2
## 3 A     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A8 CD68                 5.64               2.26                       12.8
## 4 M     M_Monocyte_CD14_resting L_NK_CD56._CD16. S100A8 ITGB2                5.73               2.14                       12.3
## 5 M     M_Monocyte_CD14_resting L_NK_CD56._CD16. S100A9 ITGB2                5.61               2.14                       12.0
## 6 A     M_Monocyte_CD14_resting L_NK_CD56._CD16. S100A8 ITGB2                5.64               2.07                       11.7
abundance_expression_info$sender_receiver_info$frq_df_group %>% head()
## # A tibble: 6 x 8
## # Groups:   group, sender [3]
##   group sender                  receiver          ligand receptor fraction_ligand_group fraction_receptor_group ligand_receptor_fraction_prod_group
##   <chr> <chr>                   <chr>             <chr>  <chr>                    <dbl>                   <dbl>                               <dbl>
## 1 A     M_Monocyte_CD16         M_Monocyte_CD16   LGALS1 PTPRC                    0.996                   0.968                               0.964
## 2 M     M_Monocyte_CD14_resting M_Monocyte_CD16   S100A8 CD68                     0.998                   0.960                               0.958
## 3 M     M_Monocyte_CD16         M_Monocyte_CD16   LGALS1 PTPRC                    0.988                   0.969                               0.957
## 4 M     M_Monocyte_CD14_resting M_Monocyte_CD16   S100A9 CD68                     0.996                   0.960                               0.957
## 5 M     M_Monocyte_CD14_resting L_T_Proliferating S100A8 ITGB2                    0.998                   0.958                               0.956
## 6 M     M_Monocyte_CD14_resting L_T_Proliferating S100A9 ITGB2                    0.996                   0.958                               0.955
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()
## # A tibble: 6 x 8
## # Groups:   group, sender [2]
##   group sender                  receiver         ligand receptor pb_ligand_group pb_receptor_group ligand_receptor_pb_prod_group
##   <chr> <fct>                   <fct>            <chr>  <chr>              <dbl>             <dbl>                         <dbl>
## 1 M     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A9 CD68                15.0             10.1                           151.
## 2 M     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A8 CD68                15.0             10.1                           150.
## 3 M     M_Monocyte_CD14_resting L_NK_CD56._CD16. S100A9 ITGB2               15.0              9.89                          148.
## 4 M     M_Monocyte_CD14_resting L_NK_CD56._CD16. S100A8 ITGB2               15.0              9.89                          148.
## 5 A     M_Monocyte_CD14_resting M_Monocyte_CD16  S100A8 CD68                14.9              9.90                          147.
## 6 A     M_Monocyte_CD14_resting L_NK_CD56._CD16. S100A8 ITGB2               14.9              9.87                          147.
```

# Step 2: Perform genome-wide differential expression analysis of receiver and sender cell types to define DE genes between the conditions of interest. Based on this analysis, we can define the logFC/p-value of ligands in senders and receptors in receivers, and define the set of affected target genes in the receiver.

Now we will go over to the multi-group, multi-sample differential
expression (DE) analysis (also called ‘differential state’ analysis by
the developers of Muscat).

### Define the contrasts and covariates of interest for the DE analysis.

Here, we want to compare each patient group to the other groups, so the
MIS-C (M) group vs healthy control siblins (S) and adult COVID19
patients (A) etcetera. To do this comparison, we need to set the
following contrasts:

``` r
contrasts_oi = c("'M-(S+A)/2','S-(M+A)/2','A-(S+M)/2'")
```

**Very Important** Note the format to indicate the contrasts! This
formatting should be adhered to very strictly, and white spaces are not
allowed! Check `?get_DE_info` for explanation about how to define this
well. The most important things are that: each contrast is surrounded by
single quotation marks, contrasts are separated by a comma without any
whitespace, and alle contrasts together are surrounded by double
quotation marks. If you compare against two groups, you should divide by
2, if you compare against three groups, you should divide by 3 etcetera.

For downstream visualizations and linking contrasts to their main group,
you need to run the following:

``` r
contrast_tbl = tibble(contrast = 
                        c("M-(S+A)/2","S-(M+A)/2", "A-(S+M)/2"), 
                      group = c("M","S","A"))
```

If you want to compare only two groups (eg M vs S), you can do like
this: `contrasts_oi = c("'M-S','S-M'")`
`contrast_tbl = tibble(contrast = c("M-S","S-M"), group = c("M","S"))`

### Perform the DE analysis for each cell type.

``` r
DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)
```

### Check DE results

Table with logFC and p-values for each gene-celltype-contrast:

``` r
DE_info$celltype_de$de_output_tidy %>% arrange(p_adj) %>% head()
## # A tibble: 6 x 9
##   gene     cluster_id              logFC logCPM     F    p_val     p_adj.loc        p_adj contrast 
##   <chr>    <chr>                   <dbl>  <dbl> <dbl>    <dbl>         <dbl>        <dbl> <chr>    
## 1 IGLV3.25 L_NK_CD56._CD16.         7.54   5.01   214 1.04e-12 0.00000000795 0.0000000383 A-(S+M)/2
## 2 TRBV11.2 L_T_TIM3._CD38._HLADR.   4.72   9.85   234 1.56e-12 0.00000000828 0.0000000573 M-(S+A)/2
## 3 IGHV3.11 L_NK_CD56._CD16.         5.72   3.09   185 8.5 e-11 0.000000325   0.00000157   A-(S+M)/2
## 4 ZDHHC1   M_Monocyte_CD14_resting  2.91   5.24   189 4.69e-11 0.000000467   0.00000173   S-(M+A)/2
## 5 IGHV1.46 L_NK_CD56._CD16.         4.58   4.13   123 2.14e-10 0.000000545   0.00000262   A-(S+M)/2
## 6 IGLV6.57 L_NK_CD56._CD16.         4.76   3.84   111 5.54e-10 0.00000106    0.0000051    A-(S+M)/2
```

It is always a good idea to check distribution of the p-values resulting
from this DE expression analysis. In order to trust the p-values, the
p-value distributions should be uniform distributions, with a peak
allowed between 0 and 0.05 if there would be a clear biological effect
in the data.

You can look at these p-value histograms in the following way:

``` r
DE_info$hist_pvals
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

The p-value distributions look OK. If this would not be the cause, this
might point to issues in the DE model definition. For example in case we
did not add all important covariates, there is substructure present in
the groups, etc.

If there are some issues, you can use the empiricall null procedure.
This is a procedure that will define empirical p-values based on the
observed distribution of the test statistic (here: logFC) and not based
on the theoretical distribution. This approach has also been used in the
Saturn package: <https://github.com/statOmics/satuRn>. We only recommend
this if the p-value distributions point to possible issues, but will
illustrate here how to do this as example:

### Empirical Null procedure

``` r
empirical_pval = TRUE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
} 
```

Table with logFC and p-values for each gene-celltype-contrast:

``` r
DE_info_emp$de_output_tidy_emp %>% arrange(p_emp) %>% head()
## # A tibble: 6 x 11
##   gene     cluster_id             logFC logCPM     F    p_val     p_adj.loc        p_adj contrast     p_emp     p_adj_emp
##   <chr>    <chr>                  <dbl>  <dbl> <dbl>    <dbl>         <dbl>        <dbl> <chr>        <dbl>         <dbl>
## 1 IGLV3.25 L_NK_CD56._CD16.        7.54   5.01   214 1.04e-12 0.00000000795 0.0000000383 A-(S+M)/2 9.90e-13 0.00000000758
## 2 TRBV11.2 L_T_TIM3._CD38._HLADR.  4.72   9.85   234 1.56e-12 0.00000000828 0.0000000573 M-(S+A)/2 1.35e-12 0.00000000720
## 3 IGHV3.11 L_NK_CD56._CD16.        5.72   3.09   185 8.5 e-11 0.000000325   0.00000157   A-(S+M)/2 8.01e-11 0.000000307  
## 4 IGHV1.46 L_NK_CD56._CD16.        4.58   4.13   123 2.14e-10 0.000000545   0.00000262   A-(S+M)/2 2.01e-10 0.000000514  
## 5 IFI27    L_T_TIM3._CD38._HLADR.  6.55   7.05    83 1.43e- 8 0.000076      0.0000389    A-(S+M)/2 2.96e-10 0.00000157   
## 6 IGLV6.57 L_NK_CD56._CD16.        4.76   3.84   111 5.54e-10 0.00000106    0.0000051    A-(S+M)/2 5.20e-10 0.000000995
```

The following plot shows the distributions of those corrected, empirical
p-values:

``` r
DE_info_emp$hist_pvals_emp
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->
The following plots show how well the correction worked. The green
fitted curve should fit well with the histogram. If not, this might
point to some issues in the DE model definition.

``` r
DE_info_emp$z_distr_plots_emp_pval
## $`L_NK_CD56._CD16..M-(S+A)/2`
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

    ## 
    ## $`L_NK_CD56._CD16..S-(M+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

    ## 
    ## $`L_NK_CD56._CD16..A-(S+M)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-3.png)<!-- -->

    ## 
    ## $`L_NK_CD56hi.M-(S+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-4.png)<!-- -->

    ## 
    ## $`L_NK_CD56hi.S-(M+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-5.png)<!-- -->

    ## 
    ## $`L_NK_CD56hi.A-(S+M)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-6.png)<!-- -->

    ## 
    ## $`L_T_Proliferating.M-(S+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-7.png)<!-- -->

    ## 
    ## $`L_T_Proliferating.S-(M+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-8.png)<!-- -->

    ## 
    ## $`L_T_Proliferating.A-(S+M)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-9.png)<!-- -->

    ## 
    ## $`L_T_TIM3._CD38._HLADR..M-(S+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-10.png)<!-- -->

    ## 
    ## $`L_T_TIM3._CD38._HLADR..S-(M+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-11.png)<!-- -->

    ## 
    ## $`L_T_TIM3._CD38._HLADR..A-(S+M)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-12.png)<!-- -->

    ## 
    ## $`M_Monocyte_CD14_resting.M-(S+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-13.png)<!-- -->

    ## 
    ## $`M_Monocyte_CD14_resting.S-(M+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-14.png)<!-- -->

    ## 
    ## $`M_Monocyte_CD14_resting.A-(S+M)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-15.png)<!-- -->

    ## 
    ## $`M_Monocyte_CD16.M-(S+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-16.png)<!-- -->

    ## 
    ## $`M_Monocyte_CD16.S-(M+A)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-17.png)<!-- -->

    ## 
    ## $`M_Monocyte_CD16.A-(S+M)/2`

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-26-18.png)<!-- -->

In general, these plots look OK.

Next, we will compare the empirical p-values to the original ones
before. We will look for the concordance between p-values rankings of
the original and empirical DE analysis (via ranking-line and upset
plots):

``` r
comparison_plots = compare_normal_emp_pvals(DE_info, DE_info_emp, adj_pval = FALSE)
comparison_plots
## [[1]]
## [[1]][[1]]
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

    ## 
    ## [[1]][[2]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->

    ## 
    ## 
    ## [[2]]
    ## [[2]][[1]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-3.png)<!-- -->

    ## 
    ## [[2]][[2]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-4.png)<!-- -->

    ## 
    ## 
    ## [[3]]
    ## [[3]][[1]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-5.png)<!-- -->

    ## 
    ## [[3]][[2]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-6.png)<!-- -->

    ## 
    ## 
    ## [[4]]
    ## [[4]][[1]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-7.png)<!-- -->

    ## 
    ## [[4]][[2]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-8.png)<!-- -->

    ## 
    ## 
    ## [[5]]
    ## [[5]][[1]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-9.png)<!-- -->

    ## 
    ## [[5]][[2]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-10.png)<!-- -->

    ## 
    ## 
    ## [[6]]
    ## [[6]][[1]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-11.png)<!-- -->

    ## 
    ## [[6]][[2]]

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-27-12.png)<!-- -->

You can see that for some cell types, the genes considered to be DE can
change.

### Conclusion of the diagnostic plots concerning the DE analysis

P-value histograms of both the normal and empirical p-values look good
here. In case the normal p-value distributions don’t look good, we
recommend to continue based on the empirical p-values.

``` r
empirical_pval = TRUE
if(empirical_pval == FALSE){
  celltype_de = DE_info$celltype_de$de_output_tidy
} else {
  celltype_de = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
}
```

In the next step, we will combine the DE information of senders and
receivers by linking their ligands and receptors together:

### Combine DE information for ligand-senders and receptors-receivers (similar to step1 - `abundance_expression_info$sender_receiver_info`)

``` r
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
```

``` r
sender_receiver_de %>% head(20)
## # A tibble: 20 x 12
##    contrast  sender                  receiver                ligand receptor lfc_ligand lfc_receptor ligand_receptor_lfc_avg p_val_ligand p_adj_ligand p_val_receptor p_adj_receptor
##    <chr>     <chr>                   <chr>                   <chr>  <chr>         <dbl>        <dbl>                   <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
##  1 M-(S+A)/2 M_Monocyte_CD16         M_Monocyte_CD16         C3     VSIG4         2.99        5.74                      4.36     0.296           1.00         0.0521           1.00 
##  2 M-(S+A)/2 M_Monocyte_CD14_resting M_Monocyte_CD16         HP     CD163         2.3         4.51                      3.40     0.00750         0.525        0.0787           1.00 
##  3 M-(S+A)/2 M_Monocyte_CD16         M_Monocyte_CD14_resting C3     VSIG4         2.99        3.64                      3.32     0.296           1.00         0.00293          0.454
##  4 M-(S+A)/2 M_Monocyte_CD14_resting M_Monocyte_CD16         C1QB   CD33          6.51        0.0355                    3.27     0.000361        0.239        0.968            1.00 
##  5 M-(S+A)/2 M_Monocyte_CD14_resting M_Monocyte_CD14_resting C1QB   CD33          6.51        0.0306                    3.27     0.000361        0.239        0.946            0.999
##  6 M-(S+A)/2 M_Monocyte_CD14_resting M_Monocyte_CD14_resting C1QB   LRP1          6.51       -0.0487                    3.23     0.000361        0.239        0.701            0.999
##  7 M-(S+A)/2 M_Monocyte_CD14_resting M_Monocyte_CD16         C1QB   LRP1          6.51       -0.316                     3.10     0.000361        0.239        0.225            1.00 
##  8 A-(S+M)/2 L_T_Proliferating       L_T_Proliferating       NCAM1  NCAM1         2.77        2.77                      2.77     0.00514         0.667        0.00514          0.667
##  9 M-(S+A)/2 M_Monocyte_CD16         M_Monocyte_CD14_resting IL1B   IL1R2         1.5         3.95                      2.72     0.136           1.00         0.00223          0.435
## 10 M-(S+A)/2 M_Monocyte_CD16         M_Monocyte_CD14_resting THBS1  SDC4          4.52        0.825                     2.67     0.109           1.00         0.126            0.888
## 11 M-(S+A)/2 M_Monocyte_CD16         M_Monocyte_CD14_resting IL1RN  IL1R2         1.24        3.95                      2.60     0.337           1.00         0.00223          0.435
## 12 M-(S+A)/2 M_Monocyte_CD16         M_Monocyte_CD16         THBS1  CD47          4.52        0.401                     2.46     0.109           1.00         0.296            1.00 
## 13 M-(S+A)/2 M_Monocyte_CD16         L_NK_CD56._CD16.        THBS1  ITGA6         4.52        0.36                      2.44     0.109           1.00         0.114            0.762
## 14 M-(S+A)/2 M_Monocyte_CD14_resting M_Monocyte_CD14_resting C1QA   CR1           4.1         0.758                     2.43     0.00166         0.394        0.0461           0.763
## 15 M-(S+A)/2 M_Monocyte_CD16         M_Monocyte_CD14_resting THBS1  CD36          4.52        0.305                     2.41     0.109           1.00         0.250            0.935
## 16 M-(S+A)/2 M_Monocyte_CD16         M_Monocyte_CD16         THBS1  ITGA6         4.52        0.303                     2.41     0.109           1.00         0.780            1.00 
## 17 M-(S+A)/2 M_Monocyte_CD16         L_NK_CD56._CD16.        THBS1  ITGA4         4.52        0.291                     2.41     0.109           1.00         0.0325           0.634
## 18 M-(S+A)/2 M_Monocyte_CD16         L_T_TIM3._CD38._HLADR.  THBS1  ITGA6         4.52        0.23                      2.38     0.109           1.00         0.426            0.900
## 19 M-(S+A)/2 M_Monocyte_CD14_resting M_Monocyte_CD14_resting IL1B   IL1R2         0.777       3.95                      2.36     0.0804          0.840        0.00223          0.435
## 20 M-(S+A)/2 M_Monocyte_CD16         M_Monocyte_CD14_resting THBS1  ITGA6         4.52        0.194                     2.36     0.109           1.00         0.613            0.997
```

# Step 3: Predict NicheNet ligand activities and NicheNet ligand-target links based on these differential expression results

## Define the parameters for the NicheNet ligand activity analysis

Here, we need to define the thresholds that will be used to consider
genes as differentially expressed or not (logFC, p-value, decision
whether to use adjusted or normal p-value, minimum fraction of cells
that should express a gene in at least one sample in a group, whether to
use the normal p-values or empirical p-values).

NicheNet ligand activity will then be calculated as the enrichment of
predicted target genes of ligands in this set of DE genes compared to
the genomic background. Here we choose for a minimum logFC of 0.50,
maximum p-value of 0.05, and minimum fraction of expression of 0.05.

``` r
logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05
```

We will here choose for applying the p-value cutoff on the normal
p-values, and not on the p-values corrected for multiple testing. This
choice was made here because this dataset has only a few samples per
group and we might have a lack of statistical power due to
pseudobulking. In case of more samples per group, and a sufficient high
number of DE genes per group-celltype (&gt; 50), we would recommend
using the adjusted p-values.

``` r
# p_val_adj = TRUE 
p_val_adj = FALSE 
```

For the NicheNet ligand-target inference, we also need to select which
top n of the predicted target genes will be considered (here: top 250
targets per ligand).

``` r
top_n_target = 250
```

The NicheNet ligand activity analysis can be run in parallel for each
receiver cell type, by changing the number of cores as defined here.
This is only recommended if you have many receiver cell type.

``` r
verbose = TRUE
cores_system = 8
n.cores = min(cores_system, sender_receiver_de$receiver %>% unique() %>% length()) # use one core per receiver cell type
```

## Run the NicheNet ligand activity analysis

(this might take some time)

``` r
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose, 
  n.cores = n.cores
)))
```

Check the DE genes used for the activity analysis

``` r
ligand_activities_targets_DEgenes$de_genes_df %>% head(20)
## # A tibble: 20 x 6
##    gene         receiver         logFC    p_val p_adj contrast 
##    <chr>        <chr>            <dbl>    <dbl> <dbl> <chr>    
##  1 TNFRSF14.AS1 L_NK_CD56._CD16. 0.957 0.0206   0.550 M-(S+A)/2
##  2 ACOT7        L_NK_CD56._CD16. 0.714 0.0113   0.493 M-(S+A)/2
##  3 KDM1A        L_NK_CD56._CD16. 0.636 0.0498   0.694 M-(S+A)/2
##  4 SMAP2        L_NK_CD56._CD16. 0.565 0.00845  0.452 M-(S+A)/2
##  5 P3H1         L_NK_CD56._CD16. 0.583 0.0196   0.544 M-(S+A)/2
##  6 CPT2         L_NK_CD56._CD16. 0.636 0.0452   0.678 M-(S+A)/2
##  7 DHCR24       L_NK_CD56._CD16. 1.15  0.00531  0.440 M-(S+A)/2
##  8 GADD45A      L_NK_CD56._CD16. 0.542 0.0414   0.663 M-(S+A)/2
##  9 SELL         L_NK_CD56._CD16. 0.583 0.00396  0.428 M-(S+A)/2
## 10 IER5         L_NK_CD56._CD16. 0.741 0.0302   0.615 M-(S+A)/2
## 11 MIR181A1HG   L_NK_CD56._CD16. 0.889 0.0307   0.622 M-(S+A)/2
## 12 CD55         L_NK_CD56._CD16. 0.584 0.0458   0.681 M-(S+A)/2
## 13 H2AW         L_NK_CD56._CD16. 0.681 0.00148  0.379 M-(S+A)/2
## 14 RNF187       L_NK_CD56._CD16. 0.666 0.00804  0.452 M-(S+A)/2
## 15 LPIN1        L_NK_CD56._CD16. 0.559 0.000491 0.342 M-(S+A)/2
## 16 POMC         L_NK_CD56._CD16. 0.595 0.0294   0.608 M-(S+A)/2
## 17 TOGARAM2     L_NK_CD56._CD16. 0.979 0.00399  0.428 M-(S+A)/2
## 18 VAMP5        L_NK_CD56._CD16. 0.596 0.000959 0.367 M-(S+A)/2
## 19 IL18R1       L_NK_CD56._CD16. 0.677 0.0494   0.694 M-(S+A)/2
## 20 CLASP1       L_NK_CD56._CD16. 0.88  0.00456  0.428 M-(S+A)/2
```

Check the output of the activity analysis

``` r
ligand_activities_targets_DEgenes$ligand_activities %>% head(20)
## # A tibble: 20 x 8
## # Groups:   receiver, contrast [1]
##    ligand activity contrast  target   ligand_target_weight receiver         direction_regulation activity_scaled
##    <chr>     <dbl> <chr>     <chr>                   <dbl> <chr>            <fct>                          <dbl>
##  1 A2M     0.00954 M-(S+A)/2 ACOT7                 0.00400 L_NK_CD56._CD16. up                             0.452
##  2 A2M     0.00954 M-(S+A)/2 CD55                  0.00377 L_NK_CD56._CD16. up                             0.452
##  3 A2M     0.00954 M-(S+A)/2 CD74                  0.00382 L_NK_CD56._CD16. up                             0.452
##  4 A2M     0.00954 M-(S+A)/2 FKBP5                 0.00447 L_NK_CD56._CD16. up                             0.452
##  5 A2M     0.00954 M-(S+A)/2 FOS                   0.00646 L_NK_CD56._CD16. up                             0.452
##  6 A2M     0.00954 M-(S+A)/2 GADD45A               0.00558 L_NK_CD56._CD16. up                             0.452
##  7 A2M     0.00954 M-(S+A)/2 ISG20                 0.00392 L_NK_CD56._CD16. up                             0.452
##  8 A2M     0.00954 M-(S+A)/2 LMNB1                 0.00383 L_NK_CD56._CD16. up                             0.452
##  9 A2M     0.00954 M-(S+A)/2 NFKB1                 0.00497 L_NK_CD56._CD16. up                             0.452
## 10 A2M     0.00954 M-(S+A)/2 NFKBIZ                0.00383 L_NK_CD56._CD16. up                             0.452
## 11 A2M     0.00954 M-(S+A)/2 PPP1R15A              0.00408 L_NK_CD56._CD16. up                             0.452
## 12 A2M     0.00954 M-(S+A)/2 SLC1A5                0.00423 L_NK_CD56._CD16. up                             0.452
## 13 A2M     0.00954 M-(S+A)/2 SMAD3                 0.00472 L_NK_CD56._CD16. up                             0.452
## 14 A2M     0.00954 M-(S+A)/2 SOCS3                 0.00488 L_NK_CD56._CD16. up                             0.452
## 15 A2M     0.00954 M-(S+A)/2 WARS1                 0.00412 L_NK_CD56._CD16. up                             0.452
## 16 AANAT   0.00319 M-(S+A)/2 FKBP5                 0.00363 L_NK_CD56._CD16. up                            -1.08 
## 17 AANAT   0.00319 M-(S+A)/2 FOS                   0.00424 L_NK_CD56._CD16. up                            -1.08 
## 18 AANAT   0.00319 M-(S+A)/2 GADD45A               0.00361 L_NK_CD56._CD16. up                            -1.08 
## 19 AANAT   0.00319 M-(S+A)/2 LMNB1                 0.00324 L_NK_CD56._CD16. up                            -1.08 
## 20 AANAT   0.00319 M-(S+A)/2 SLC16A3               0.00311 L_NK_CD56._CD16. up                            -1.08
```

# Step 4: Use the information collected above to prioritize all sender-ligand—receiver-receptor pairs.

In the 3 previous steps, we calculated expression, differential
expression and NicheNet activity information. Now we will combine these
different types of information in one prioritization scheme.

MultiNicheNet allows the user to define the weights of the following
criteria to prioritize ligand-receptor interactions:

-   Upregulation of the ligand in a sender cell type and/or upregulation
    of the receptor in a receiver cell type - in the condition of
    interest. : `scaled_lfc_ligand`, `scaled_p_val_ligand`,
    `scaled_lfc_receptor`, and `scaled_p_val_receptor`
-   Sufficiently high expression levels of ligand and receptor in many
    samples of the same group (to mitigate the influence of outlier
    samples). : `fraction_expressing_ligand_receptor`
-   Cell-type and condition specific expression of the ligand in the
    sender cell type and receptor in the receiver cell type (to mitigate
    the influence of upregulated but still relatively weakly expressed
    ligands/receptors) : `scaled_avg_exprs_ligand`,
    `scaled_avg_frq_ligand`, `scaled_avg_exprs_receptor`, and
    `scaled_avg_frq_receptor`
-   High NicheNet ligand activity, to further prioritize ligand-receptor
    pairs based on their predicted effect of the ligand-receptor
    interaction on the gene expression in the receiver cell type :
    `scaled_activity` and `scaled_activity_scaled`
    (scaled\_activity=absolute value of ligand activity;
    scaled\_activity\_scaled=scaled ligand activity value that is
    comparable between different receiver settings - recommended to put
    more weight on this)
-   High relative abundance of sender and/or receiver in the condition
    of interest: `scaled_abundance_sender` and
    `scaled_abundance_receiver`

The different properties of the sender-ligand—receiver-receptor pairs
can be weighted according to the user’s preference and insight in the
dataset at hand.

## Define the prioritization weights, and prepare grouping objects

We will set our preference for this dataset as follows - and recommend
the user to use the same weights by default:

``` r
prioritizing_weights_DE = c("de_ligand" = 1,
                         "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                         "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                         "abund_receiver" = 0)
```

``` r
prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)
```

Make necessary grouping data frame

``` r
sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}
```

Crucial note: grouping\_tbl: group should be the same as in the
contrast\_tbl, and as in the expression info tables! Rename accordingly
if this would not be the case. If you followed the guidelines of this
tutorial closely, there should be no problem.

## Run the prioritization

``` r
prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))
```

Check the output tables

First: group-based summary table

``` r
prioritization_tables$group_prioritization_tbl %>% head(20)
## # A tibble: 20 x 60
##    contrast  group sender  receiver  ligand  receptor lfc_ligand lfc_receptor ligand_receptor_~ p_val_ligand p_adj_ligand p_val_receptor p_adj_receptor activity direction_regula~ activity_scaled lr_interaction id      avg_ligand_group avg_receptor_gr~ ligand_receptor~ fraction_ligand~ fraction_recept~ ligand_receptor~ rel_abundance_s~ rel_abundance_s~ sender_receiver~
##    <chr>     <chr> <chr>   <chr>     <chr>   <chr>         <dbl>        <dbl>             <dbl>        <dbl>        <dbl>          <dbl>          <dbl>    <dbl> <fct>                       <dbl> <chr>          <chr>              <dbl>            <dbl>            <dbl>            <dbl>            <dbl>            <dbl>            <dbl>            <dbl>            <dbl>
##  1 A-(S+M)/2 A     M_Mono~ L_T_TIM3~ NAMPT   ITGB1         1.25         0.831             1.04     0.0000663       0.0412       0.00374           0.273 0.00308  up                        -0.170  NAMPT_ITGB1    NAMPT_~            1.11             0.907           1.01              0.558            0.554           0.309            0.212             0.221            0.217
##  2 A-(S+M)/2 A     M_Mono~ L_T_TIM3~ NAMPT   ITGB1         1.25         0.831             1.04     0.0000663       0.0412       0.00374           0.273 0.0188   down                       3.39   NAMPT_ITGB1    NAMPT_~            1.11             0.907           1.01              0.558            0.554           0.309            0.212             0.221            0.217
##  3 A-(S+M)/2 A     M_Mono~ L_T_TIM3~ NAMPT   ITGB1         0.622        0.831             0.726    0.0991          0.949        0.00374           0.273 0.00308  up                        -0.170  NAMPT_ITGB1    NAMPT_~            1.35             0.907           1.22              0.730            0.554           0.405            0.0441            0.221            0.133
##  4 A-(S+M)/2 A     M_Mono~ L_T_TIM3~ NAMPT   ITGB1         0.622        0.831             0.726    0.0991          0.949        0.00374           0.273 0.0188   down                       3.39   NAMPT_ITGB1    NAMPT_~            1.35             0.907           1.22              0.730            0.554           0.405            0.0441            0.221            0.133
##  5 M-(S+A)/2 M     L_T_Pr~ M_Monocy~ IFNG    IFNGR2        1.22         0.389             0.804    0.00699         0.481        0.104             0.881 0.0199   up                         2.86   IFNG_IFNGR2    IFNG_I~            0.149            0.931           0.139             0.152            0.547           0.0829           0.945             0.358            0.652
##  6 M-(S+A)/2 M     L_T_Pr~ M_Monocy~ IFNG    IFNGR2        1.22         0.389             0.804    0.00699         0.481        0.104             0.881 0.000600 down                       0.294  IFNG_IFNGR2    IFNG_I~            0.149            0.931           0.139             0.152            0.547           0.0829           0.945             0.358            0.652
##  7 M-(S+A)/2 M     L_T_Pr~ M_Monocy~ IFNG    IFNGR1        1.22         0.515             0.867    0.00699         0.481        0.281             1.00  0.0183   up                         2.91   IFNG_IFNGR1    IFNG_I~            0.149            0.575           0.0856            0.152            0.498           0.0755           0.945             0.442            0.693
##  8 M-(S+A)/2 M     L_T_Pr~ M_Monocy~ IFNG    IFNGR1        1.22         0.515             0.867    0.00699         0.481        0.281             1.00  0.00345  down                      -1.87   IFNG_IFNGR1    IFNG_I~            0.149            0.575           0.0856            0.152            0.498           0.0755           0.945             0.442            0.693
##  9 A-(S+M)/2 A     M_Mono~ L_T_Prol~ NAMPT   ITGB1         1.25         0.551             0.901    0.0000663       0.0412       0.0528            0.810 0.00881  up                         0.0437 NAMPT_ITGB1    NAMPT_~            1.11             0.807           0.895             0.558            0.552           0.308            0.212             0.270            0.241
## 10 A-(S+M)/2 A     M_Mono~ L_T_Prol~ NAMPT   ITGB1         1.25         0.551             0.901    0.0000663       0.0412       0.0528            0.810 0.00620  down                       3.14   NAMPT_ITGB1    NAMPT_~            1.11             0.807           0.895             0.558            0.552           0.308            0.212             0.270            0.241
## 11 M-(S+A)/2 M     L_T_TI~ M_Monocy~ IFNG    IFNGR2        1.25         0.389             0.820    0.0310          0.364        0.104             0.881 0.0199   up                         2.86   IFNG_IFNGR2    IFNG_I~            0.155            0.931           0.144             0.114            0.547           0.0623           1.00              0.358            0.679
## 12 M-(S+A)/2 M     L_T_TI~ M_Monocy~ IFNG    IFNGR2        1.25         0.389             0.820    0.0310          0.364        0.104             0.881 0.000600 down                       0.294  IFNG_IFNGR2    IFNG_I~            0.155            0.931           0.144             0.114            0.547           0.0623           1.00              0.358            0.679
## 13 M-(S+A)/2 M     L_T_TI~ M_Monocy~ IFNG    IFNGR1        1.25         0.515             0.883    0.0310          0.364        0.281             1.00  0.0183   up                         2.91   IFNG_IFNGR1    IFNG_I~            0.155            0.575           0.0889            0.114            0.498           0.0567           1.00              0.442            0.721
## 14 M-(S+A)/2 M     L_T_TI~ M_Monocy~ IFNG    IFNGR1        1.25         0.515             0.883    0.0310          0.364        0.281             1.00  0.00345  down                      -1.87   IFNG_IFNGR1    IFNG_I~            0.155            0.575           0.0889            0.114            0.498           0.0567           1.00              0.442            0.721
## 15 S-(M+A)/2 S     M_Mono~ M_Monocy~ HLA.DMA CD4           0.88         0.859             0.869    0.000153        0.0950       0.00678           0.381 0.0137   up                         1.80   HLA.DMA_CD4    HLA.DM~            0.760            0.520           0.396             0.442            0.330           0.146            0.653             0.653            0.653
## 16 S-(M+A)/2 S     M_Mono~ M_Monocy~ HLA.DMA CD4           0.88         0.859             0.869    0.000153        0.0950       0.00678           0.381 0.0137   down                      -0.621  HLA.DMA_CD4    HLA.DM~            0.760            0.520           0.396             0.442            0.330           0.146            0.653             0.653            0.653
## 17 M-(S+A)/2 M     M_Mono~ L_T_Prol~ HLA.DRA LAG3          0.598        1.75              1.17     0.0766          1.00         0.000113          0.152 0.0205   up                         2.68   HLA.DRA_LAG3   HLA.DR~            2.92             0.669           1.95              0.990            0.520           0.515            0.442             0.945            0.693
## 18 M-(S+A)/2 M     M_Mono~ L_T_Prol~ HLA.DRA LAG3          0.598        1.75              1.17     0.0766          1.00         0.000113          0.152 0.00372  down                      -0.147  HLA.DRA_LAG3   HLA.DR~            2.92             0.669           1.95              0.990            0.520           0.515            0.442             0.945            0.693
## 19 M-(S+A)/2 M     L_T_Pr~ M_Monocy~ IFNG    IFNGR1        1.22         0.265             0.742    0.00699         0.481        0.355             0.960 0.0199   up                         2.86   IFNG_IFNGR1    IFNG_I~            0.149            0.424           0.0632            0.152            0.293           0.0444           0.945             0.358            0.652
## 20 M-(S+A)/2 M     L_T_Pr~ M_Monocy~ IFNG    IFNGR1        1.22         0.265             0.742    0.00699         0.481        0.355             0.960 0.000600 down                       0.294  IFNG_IFNGR1    IFNG_I~            0.149            0.424           0.0632            0.152            0.293           0.0444           0.945             0.358            0.652
## # ... with 33 more variables: lfc_pval_ligand <dbl>, p_val_ligand_adapted <dbl>, scaled_lfc_ligand <dbl>, scaled_p_val_ligand <dbl>, scaled_lfc_pval_ligand <dbl>, scaled_p_val_ligand_adapted <dbl>, lfc_pval_receptor <dbl>, p_val_receptor_adapted <dbl>, scaled_lfc_receptor <dbl>, scaled_p_val_receptor <dbl>, scaled_lfc_pval_receptor <dbl>,
## #   scaled_p_val_receptor_adapted <dbl>, activity_up <dbl>, activity_scaled_up <dbl>, scaled_activity_scaled_up <dbl>, scaled_activity_up <dbl>, activity_down <dbl>, activity_scaled_down <dbl>, scaled_activity_scaled_down <dbl>, scaled_activity_down <dbl>, scaled_avg_exprs_ligand <dbl>, scaled_avg_frq_ligand <dbl>, pb_ligand_group <dbl>, scaled_pb_ligand <dbl>,
## #   scaled_avg_exprs_receptor <dbl>, scaled_avg_frq_receptor <dbl>, pb_receptor_group <dbl>, scaled_pb_receptor <dbl>, fraction_expressing_ligand_receptor <dbl>, max_scaled_activity <dbl>, na.rm <lgl>, prioritization_score <dbl>, top_group <chr>
```

Second: sample-based summary table: contains expression information of
each LR pair per sample

``` r
prioritization_tables$sample_prioritization_tbl %>% head(20)
## # A tibble: 20 x 26
##    sample sender                  receiver                ligand receptor avg_ligand avg_receptor ligand_receptor_~ fraction_ligand fraction_recept~ ligand_receptor~ pb_ligand pb_receptor ligand_receptor~ group prioritization_~ lr_interaction id    scaled_LR_prod scaled_LR_frac scaled_LR_pb_pr~ n_cells_receiver keep_receiver n_cells_sender keep_sender keep_sender_rec~
##    <chr>  <chr>                   <chr>                   <chr>  <chr>         <dbl>        <dbl>             <dbl>           <dbl>            <dbl>            <dbl>     <dbl>       <dbl>            <dbl> <chr>            <dbl> <chr>          <chr>          <dbl>          <dbl>            <dbl>            <dbl>         <dbl>          <dbl>       <dbl> <fct>           
##  1 M4     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A8 CD68           6.48         2.81              18.2           1                1                1          16.1       11.1              179. M                1.18  S100A8_CD68    S100~          1.94           0.710            1.95               224             1            379           1 Sender & Receiv~
##  2 A1     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A8 CD68           6.33         2.82              17.9           0.984            1                0.984      15.8       10.3              163. A                1.12  S100A8_CD68    S100~          1.83           0.526            1.01                 6             0             61           1 Receiver absent 
##  3 M4     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A9 CD68           6.00         2.81              16.9           1                1                1          15.8       11.1              174. M                1.20  S100A9_CD68    S100~          1.98           0.729            1.84               224             1            379           1 Sender & Receiv~
##  4 M4     L_NK_CD56._CD16.        M_Monocyte_CD16         B2M    LILRB1         5.07         3.01              15.3           1                0.996            0.996      14.0       11.5              161. M                0.960 B2M_LILRB1     B2M_~          2.33           1.11             2.27               224             1            573           1 Sender & Receiv~
##  5 A1     M_Monocyte_CD14_resting L_NK_CD56._CD16.        S100A8 ITGB2          6.33         2.41              15.2           0.984            0.930            0.915      15.8       10.3              162. A                1.01  S100A8_ITGB2   S100~          1.85           0.754            1.78              1302             1             61           1 Sender & Receiv~
##  6 M5     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A8 CD68           5.85         2.56              15.0           1                0.979            0.979      15.0       10.8              162. M                1.18  S100A8_CD68    S100~          0.911          0.477            0.975               96             1             54           1 Sender & Receiv~
##  7 M4     M_Monocyte_CD16         M_Monocyte_CD16         B2M    LILRB1         4.97         3.01              15.0           1                0.996            0.996      14.2       11.5              162. M                1.10  B2M_LILRB1     B2M_~          2.43           1.10             2.30               224             1            224           1 Sender & Receiv~
##  8 M4     M_Monocyte_CD14_resting L_NK_CD56._CD16.        S100A8 ITGB2          6.48         2.30              14.9           1                0.941            0.941      16.1       10.1              164. M                1.10  S100A8_ITGB2   S100~          1.71           1.15             1.88               573             1            379           1 Sender & Receiv~
##  9 M4     L_T_TIM3._CD38._HLADR.  M_Monocyte_CD16         B2M    LILRB1         4.95         3.01              14.9           1                0.996            0.996      13.8       11.5              158. M                1.00  B2M_LILRB1     B2M_~          2.31           1.10             2.23               224             1            130           1 Sender & Receiv~
## 10 M5     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A9 CD68           5.75         2.56              14.7           1                0.979            0.979      15.0       10.8              162. M                1.20  S100A9_CD68    S100~          1.16           0.495            1.07                96             1             54           1 Sender & Receiv~
## 11 M7     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A8 CD68           5.50         2.63              14.5           0.995            0.985            0.980      14.7       10.8              159. M                1.18  S100A8_CD68    S100~          0.734          0.480            0.770              196             1            190           1 Sender & Receiv~
## 12 A1     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A9 CD68           5.10         2.82              14.4           0.984            1                0.984      14.5       10.3              149. A                0.921 S100A9_CD68    S100~          1.03           0.545            0.248                6             0             61           1 Receiver absent 
## 13 M4     L_T_Proliferating       M_Monocyte_CD16         B2M    LILRB1         4.73         3.01              14.2           1                0.996            0.996      13.5       11.5              154. M                0.929 B2M_LILRB1     B2M_~          2.24           1.10             2.15               224             1             84           1 Sender & Receiv~
## 14 M4     L_NK_CD56hi             M_Monocyte_CD16         B2M    LILRB1         4.71         3.01              14.2           0.982            0.996            0.978      13.7       11.5              157. M                0.861 B2M_LILRB1     B2M_~          2.29           0.989            2.36               224             1             56           1 Sender & Receiv~
## 15 M4     M_Monocyte_CD14_resting M_Monocyte_CD14_resting S100A8 ITGB2          6.48         2.18              14.1           1                0.947            0.947      16.1       10.1              163. M                1.18  S100A8_ITGB2   S100~          1.81           1.52             1.84               379             1            379           1 Sender & Receiv~
## 16 M4     M_Monocyte_CD14_resting L_NK_CD56._CD16.        S100A9 ITGB2          6.00         2.30              13.8           1                0.941            0.941      15.8       10.1              160. M                1.11  S100A9_ITGB2   S100~          1.63           1.16             1.93               573             1            379           1 Sender & Receiv~
## 17 M6     M_Monocyte_CD14_resting L_NK_CD56._CD16.        S100A8 ITGB2          6.48         2.11              13.7           1                0.915            0.915      16.0        9.80             157. M                1.10  S100A8_ITGB2   S100~          1.16           0.754            1.30                47             1             45           1 Sender & Receiv~
## 18 M3     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A9 CD68           5.67         2.41              13.7           1                0.8              0.8        15.0       10.3              154. M                1.20  S100A9_CD68    S100~          0.739         -1.52             0.536                5             0           1120           1 Receiver absent 
## 19 M7     M_Monocyte_CD14_resting M_Monocyte_CD16         S100A9 CD68           5.17         2.63              13.6           0.984            0.985            0.969      14.6       10.8              158. M                1.20  S100A9_CD68    S100~          0.719          0.382            0.804              196             1            190           1 Sender & Receiv~
## 20 A1     M_Monocyte_CD14_resting L_T_Proliferating       S100A8 ITGB2          6.33         2.14              13.6           0.984            0.962            0.946      15.8        9.64             152. A                0.972 S100A8_ITGB2   S100~          1.83           0.517            1.57               104             1             61           1 Sender & Receiv~
```

# Step 5: Add information on prior knowledge and expression correlation between LR and target expression.

In multi-sample datasets, we have the opportunity to look whether
expression of ligand-receptor across all samples is correlated with the
expression of their by NicheNet predicted target genes. This is what we
will do with the following line of code:

``` r
lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)
```

# Save all the output of MultiNicheNet

To avoid needing to redo the analysis later. All the output written down
here is sufficient to make all in-built downstream visualizations.

``` r
path = "./"

multinichenet_output = list(
    celltype_info = abundance_expression_info$celltype_info,
    celltype_de = celltype_de,
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  ) 
multinichenet_output = make_lite_output(multinichenet_output)

save = FALSE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output.rds"))

}
```

# Visualization of the results of the cell-cell communication analysis

In a first instance, we will look at the broad overview of prioritized
interactions via condition-specific Circos plots.

## Circos plot of top-prioritized links

We will look here at the top 50 predictions across all contrasts,
senders, and receivers of interest.

``` r
prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, rank_per_group = FALSE)
```

``` r
prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-47-2.png)<!-- -->![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-47-3.png)<!-- -->![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-47-4.png)<!-- -->

Now you can also make a full circos plot for one group of interest,
where will show the top30 per group

``` r
prioritized_tbl_oi_A_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = "A")
prioritized_tbl_oi_M_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = "M")
prioritized_tbl_oi_S_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = "S")
```

``` r
circos_A = make_circos_one_group(prioritized_tbl_oi_A_30, colors_sender, colors_receiver)
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-49-2.png)<!-- -->

``` r
circos_M = make_circos_one_group(prioritized_tbl_oi_M_30, colors_sender, colors_receiver)
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-49-3.png)<!-- -->![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-49-4.png)<!-- -->

``` r
circos_S = make_circos_one_group(prioritized_tbl_oi_S_30, colors_sender, colors_receiver)
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-49-5.png)<!-- -->![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-49-6.png)<!-- -->

## Visualization of scaled ligand-receptor pseudobulk products and ligand activity

Now we will visualize per sample the scaled product of ligand and
receptor expression. Samples that were left out of the DE analysis are
indicated with a smaller dot (this helps to indicate the samples that
did not contribute to the calculation of the logFC, and thus not
contributed to the final prioritization)

We will now check the top 50 interactions specific for the MIS-C group

``` r
group_oi = "M"
```

``` r
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, groups_oi = group_oi)

plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

Typically, there are way more than 50 differentially expressed and
active ligand-receptor pairs per group across all sender-receiver
combinations. Therefore it might be useful to zoom in on specific cell
types as senders/receivers:

Eg M\_Monocyte\_CD16 as receiver:

``` r
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, groups_oi = group_oi, receivers_oi = "M_Monocyte_CD16")

plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

Eg M\_Monocyte\_CD16 as sender:

``` r
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, groups_oi = group_oi, senders_oi = "M_Monocyte_CD16")

plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

Now the entire top50 for the S group

``` r
group_oi = "S"
```

``` r
prioritized_tbl_oi_S_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, groups_oi = group_oi)

plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_S_50)
plot_oi
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

Now the entire top50 for the A group:

``` r
group_oi = "A"
```

``` r
prioritized_tbl_oi_A_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, groups_oi = group_oi)

plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_A_50)
plot_oi
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

**Note**: Use `make_sample_lr_prod_activity_batch_plots` if you have
batches and want to visualize them on this plot!

## Visualization of ligand-activity, ligand-target links, and target gene expression

In another type of plot, we can visualize the ligand activities for a
group-receiver combination, and show the predicted ligand-target links,
and also the expression of the predicted target genes across samples.

For this, we now need to define a receiver cell type of interest. As
example, we will take `L_T_TIM3._CD38._HLADR.` cells as receiver, and
look at the top 20 senderLigand-receiverReceptor pairs with these cells
as receiver.

``` r
group_oi = "M"
receiver_oi = "L_T_TIM3._CD38._HLADR."
prioritized_tbl_oi_M_10 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 10, groups_oi = group_oi, receivers_oi = receiver_oi)
```

``` r
combined_plot = make_ligand_activity_target_plot(group_oi, receiver_oi, prioritized_tbl_oi_M_10, multinichenet_output$prioritization_tables, multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, ligand_target_matrix, plot_legend = FALSE)
combined_plot
## $combined_plot
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

    ## 
    ## $legends

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-59-2.png)<!-- -->

Of course you can look at other receivers as well:

``` r
group_oi = "M"
receiver_oi = "M_Monocyte_CD16"
prioritized_tbl_oi_M_10 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 10, groups_oi = group_oi, receivers_oi = receiver_oi)
```

``` r
combined_plot = make_ligand_activity_target_plot(group_oi, receiver_oi, prioritized_tbl_oi_M_10, multinichenet_output$prioritization_tables, multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, ligand_target_matrix, plot_legend = FALSE)
combined_plot
## $combined_plot
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

    ## 
    ## $legends

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-61-2.png)<!-- -->
\#\# Show DE genes only:

``` r
group_oi = "M"
receiver_oi = "M_Monocyte_CD16"
targets_oi = multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% arrange(p_val) %>% filter(receiver == receiver_oi & logFC > 1 & contrast == contrast_tbl %>% filter(group == group_oi) %>% pull(contrast)) %>% pull(gene) %>% unique()

p_target = make_DEgene_dotplot_pseudobulk(genes_oi = targets_oi, celltype_info = multinichenet_output$celltype_info, prioritization_tables = multinichenet_output$prioritization_tables, celltype_oi = receiver_oi, multinichenet_output$grouping_tbl)
p_target$pseudobulk_plot
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

``` r
p_target$singlecell_plot
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-62-2.png)<!-- -->

**Note** Use `make_DEgene_dotplot_pseudobulk_batch` if you want to
indicate the batch of each sample to the plot

## Visualization of expression-correlated target genes of ligand-receptor pairs

Before, we had calculated the correlation in expression between
ligand-receptor pairs and DE genes. Now we will filter out correlated
ligand-receptor –&gt; target links that both show high expression
correlation (spearman or activity &gt; 0.50 in this example) and have
some prior knowledge to support their link.

``` r
group_oi = "M"
receiver_oi = "M_Monocyte_CD16"
lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>% inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi, receiver == receiver_oi)
lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50)) # downregulation -- negative correlation
lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)
```

Now we will visualize the top correlated hits for the LR pairs that are
also in the top 50 LR pairs discriminating the groups from each other:

``` r
prioritized_tbl_oi = get_top_n_lr_pairs(prioritization_tables, 50, groups_oi = group_oi, receivers_oi = receiver_oi)
```

First: show the LR–&gt;Target heatmap of prior knowledge supported and
correlated links (shows both measure of correlation and measure of prior
knowledge)

For upregulated targets:

``` r
lr_target_prior_cor_heatmap = make_lr_target_prior_cor_heatmap(lr_target_prior_cor_filtered_up %>% filter(receiver == receiver_oi & id %in% prioritized_tbl_oi$id), add_grid = TRUE)
lr_target_prior_cor_heatmap
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

For downregulated targets:

``` r
lr_target_prior_cor_heatmap = make_lr_target_prior_cor_heatmap(lr_target_prior_cor_filtered_down %>% filter(receiver == receiver_oi & id %in% prioritized_tbl_oi$id), add_grid = TRUE)
lr_target_prior_cor_heatmap
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

Now visualize these links together with the LR and target expression, to
visualize the expression correlation

``` r
lr_target_correlation_plot = make_lr_target_correlation_plot(multinichenet_output$prioritization_tables, prioritized_tbl_oi,  lr_target_prior_cor_filtered , multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, receiver_oi,plot_legend = FALSE)
lr_target_correlation_plot$combined_plot
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

You can also visualize the expression correlation in the following way
for a selected LR pair and their targets:

``` r
ligand_oi = "IFNG"
receptor_oi = "IFNGR2"
sender_oi = "L_T_Proliferating"
receiver_oi = "M_Monocyte_CD16"
lr_target_scatter_plot = make_lr_target_scatter_plot(multinichenet_output$prioritization_tables, ligand_oi, receptor_oi, sender_oi, receiver_oi, multinichenet_output$celltype_info, multinichenet_output$grouping_tbl, lr_target_prior_cor_filtered)
lr_target_scatter_plot
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->

``` r
ligand_oi = "IFNG"
receptor_oi = "IFNGR1"
sender_oi = "L_T_TIM3._CD38._HLADR."
receiver_oi = "M_Monocyte_CD16"
lr_target_scatter_plot = make_lr_target_scatter_plot(multinichenet_output$prioritization_tables, ligand_oi, receptor_oi, sender_oi, receiver_oi, multinichenet_output$celltype_info, multinichenet_output$grouping_tbl, lr_target_prior_cor_filtered)
lr_target_scatter_plot
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

For these targets, you can also visualize the ‘prior knowledge’
ligand-receptor-to-target signaling paths. This is done similarly to the
workflow described in
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_target_signaling_path.md>

``` r
organism = "human"
if(organism == "human"){
  sig_network = readRDS(url("https://zenodo.org/record/5884439/files/signaling_network_human_21122021.rds")) %>% mutate(from = make.names(from), to = make.names(to))
  gr_network = readRDS(url("https://zenodo.org/record/5884439/files/gr_network_human_21122021.rds")) %>% mutate(from = make.names(from), to = make.names(to))
  ligand_tf_matrix = readRDS(url("https://zenodo.org/record/5884439/files/ligand_tf_matrix_nsga2r_final.rds"))
  colnames(ligand_tf_matrix) = colnames(ligand_tf_matrix) %>% make.names()
  rownames(ligand_tf_matrix) = rownames(ligand_tf_matrix) %>% make.names()
  weighted_networks = readRDS(url("https://zenodo.org/record/5884439/files/weighted_networks_nsga2r_final.rds"))
  weighted_networks$lr_sig = weighted_networks$lr_sig %>% mutate(from = make.names(from), to = make.names(to))
  weighted_networks$gr = weighted_networks$gr %>% mutate(from = make.names(from), to = make.names(to))
} else if(organism == "mouse"){
  sig_network = readRDS(url("https://zenodo.org/record/5884439/files/signaling_network_mouse_21122021.rds")) %>% mutate(from = make.names(from), to = make.names(to))
  gr_network = readRDS(url("https://zenodo.org/record/5884439/files/gr_network_mouse_21122021.rds")) %>% mutate(from = make.names(from), to = make.names(to))
  ligand_tf_matrix = readRDS(url("https://zenodo.org/record/5884439/files/ligand_tf_matrix_nsga2r_final_mouse.rds"))
  colnames(ligand_tf_matrix) = colnames(ligand_tf_matrix) %>% make.names()
  rownames(ligand_tf_matrix) = rownames(ligand_tf_matrix) %>% make.names()
  weighted_networks = readRDS(url("https://zenodo.org/record/5884439/files/weighted_networks_nsga2r_final_mouse.rds"))
  weighted_networks$lr_sig = weighted_networks$lr_sig %>% mutate(from = make.names(from), to = make.names(to))
  weighted_networks$gr = weighted_networks$gr %>% mutate(from = make.names(from), to = make.names(to))
}
```

``` r
targets_all = lr_target_prior_cor_filtered %>% filter(ligand == ligand_oi & receiver == receiver_oi & sender == sender_oi & receptor == receptor_oi)  %>% pull(target) %>% unique()
  
active_signaling_network = nichenetr::get_ligand_signaling_path_with_receptor(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligand_oi, receptors_all = receptor_oi, targets_all = targets_all, weighted_networks = weighted_networks, top_n_regulators = 2)
data_source_network = nichenetr::infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
  
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
colors = c("ligand" = "purple", "receptor" = "orange", "target" = "royalblue", "mediator" = "grey60")
ggraph_signaling_path = suppressWarnings(make_ggraph_signaling_path(active_signaling_network_min_max, colors, ligand_oi, receptor_oi, targets_all))
ggraph_signaling_path$plot
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

``` r
data_source_network %>% head()
## # A tibble: 6 x 7
##   from  to      source                        database        layer      ligand receptor
##   <chr> <chr>   <chr>                         <chr>           <chr>      <chr>  <chr>   
## 1 CD40  TNFAIP2 CytoSig_all                   CytoSig         regulatory <NA>   <NA>    
## 2 CD40  TNFAIP2 CytoSig_signature             CytoSig         regulatory <NA>   <NA>    
## 3 IFNG  B2M     lr_evex_regulation_expression evex_expression regulatory <NA>   <NA>    
## 4 IFNG  B2M     CytoSig_signature             CytoSig         regulatory <NA>   <NA>    
## 5 IFNG  B2M     NicheNet_LT_frequent          NicheNet_LT     regulatory <NA>   <NA>    
## 6 IFNG  C1QB    CytoSig_all                   CytoSig         regulatory <NA>   <NA>
```

## Intercellular regulatory network systems view

In the plots before, we demonstrated that some DE genes have both
expression correlation and prior knowledge support to be downstream of
ligand-receptor pairs. Interestingly, some target genes can be ligands
or receptors themselves. This illustrates that cells can send signals to
other cells, who as a response to these signals produce signals
themselves to feedback to the original sender cells, or who will effect
other cell types.

As last plot, we can generate a ‘systems’ view of these intercellular
feedback and cascade processes than can be occuring between the
different cell populations involved. In this plot, we will draw links
between ligands of sender cell types their ligand/receptor-annotated
target genes in receiver cell types. So links are ligand-target links (=
gene regulatory links) and not ligand-receptor protein-protein
interactions!

``` r
prioritized_tbl_oi = get_top_n_lr_pairs(prioritization_tables, 150, rank_per_group = FALSE)

lr_target_prior_cor_filtered = prioritization_tables$group_prioritization_tbl$group %>% unique() %>% lapply(function(group_oi){
  lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>% inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi)
  lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
  lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50))
  lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)
}) %>% bind_rows()
```

``` r
colors_sender["L_T_TIM3._CD38._HLADR."] = "pink" # the  original lightgreen with white font is not very readable
graph_plot = make_ggraph_ligand_target_links(lr_target_prior_cor_filtered = lr_target_prior_cor_filtered, prioritized_tbl_oi = prioritized_tbl_oi, colors = colors_sender)
graph_plot$plot
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->

``` r
graph_plot$source_df_lt %>% head()
## # A tibble: 6 x 6
##   sender                       receiver                        direction_regulation group type          weight
##   <chr>                        <chr>                           <fct>                <chr> <chr>          <dbl>
## 1 M_Monocyte_CD16_EBI3         L_T_TIM3._CD38._HLADR._LGALS3BP up                   A     Ligand-Target      1
## 2 L_NK_CD56._CD16._CCL4        M_Monocyte_CD14_resting_ADM     up                   A     Ligand-Target      1
## 3 L_T_Proliferating_CCL4       M_Monocyte_CD14_resting_ADM     up                   A     Ligand-Target      1
## 4 M_Monocyte_CD14_resting_CD55 M_Monocyte_CD14_resting_ADM     up                   A     Ligand-Target      1
## 5 L_NK_CD56hi_CD55             M_Monocyte_CD14_resting_ADM     up                   A     Ligand-Target      1
## 6 L_NK_CD56._CD16._CCL3        M_Monocyte_CD14_resting_ADM     up                   A     Ligand-Target      1
graph_plot$nodes_df %>% head()
##                                                      node                celltype gene       type_gene
## L_T_TIM3._CD38._HLADR._CD28   L_T_TIM3._CD38._HLADR._CD28  L_T_TIM3._CD38._HLADR. CD28 ligand/receptor
## M_Monocyte_CD14_resting_CD86 M_Monocyte_CD14_resting_CD86 M_Monocyte_CD14_resting CD86 ligand/receptor
## M_Monocyte_CD16_EBI3                 M_Monocyte_CD16_EBI3         M_Monocyte_CD16 EBI3          ligand
## L_NK_CD56._CD16._CCL4               L_NK_CD56._CD16._CCL4        L_NK_CD56._CD16. CCL4          ligand
## L_T_Proliferating_CCL4             L_T_Proliferating_CCL4       L_T_Proliferating CCL4          ligand
## M_Monocyte_CD14_resting_CD55 M_Monocyte_CD14_resting_CD55 M_Monocyte_CD14_resting CD55          ligand
```

## Show ligand activities for each receiver-group combination

In the next type of plot, we plot all the ligand activities (Both scaled
and absolute activities) of each receiver-group combination. This can
give us some insights in active signaling pathways across groups. Note
that we can thus show top ligands based on ligand activity - agnostic of
expression in sender.

``` r
ligands_oi = multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl %>% inner_join(contrast_tbl) %>% 
  group_by(group, receiver) %>% distinct(ligand, receiver, group, activity) %>% 
  top_n(5, activity) %>% pull(ligand) %>% unique()

plot_oi = make_ligand_activity_plots(multinichenet_output$prioritization_tables, ligands_oi, contrast_tbl, widths = NULL)
plot_oi
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

## Zoom in on specific ligand-receptor interactions: show their expression in the single-cell data!

Single-cell-based Feature, and Violin plots of ligand-receptor
interaction of interest: `make_ligand_receptor_feature_plot` and
`make_ligand_receptor_violin_plot`

It is often useful to zoom in on specific ligand-receptor interactions
of interest by looking in more detail to their expression at the single
cell level

Check the highest scoring links based on the general prioritization
score. Here we will pick one of those to visualize.

``` r
prioritized_tbl_oi %>% head(10)
## # A tibble: 10 x 8
##    group sender                  receiver                ligand  receptor id                                                          prioritization_score prioritization_rank
##    <chr> <chr>                   <chr>                   <chr>   <chr>    <chr>                                                                      <dbl>               <dbl>
##  1 A     M_Monocyte_CD14_resting L_T_TIM3._CD38._HLADR.  NAMPT   ITGB1    NAMPT_ITGB1_M_Monocyte_CD14_resting_L_T_TIM3._CD38._HLADR.                  1.36                   1
##  2 A     M_Monocyte_CD16         L_T_TIM3._CD38._HLADR.  NAMPT   ITGB1    NAMPT_ITGB1_M_Monocyte_CD16_L_T_TIM3._CD38._HLADR.                          1.34                   2
##  3 M     L_T_Proliferating       M_Monocyte_CD14_resting IFNG    IFNGR2   IFNG_IFNGR2_L_T_Proliferating_M_Monocyte_CD14_resting                       1.31                   3
##  4 M     L_T_Proliferating       M_Monocyte_CD16         IFNG    IFNGR1   IFNG_IFNGR1_L_T_Proliferating_M_Monocyte_CD16                               1.31                   4
##  5 A     M_Monocyte_CD14_resting L_T_Proliferating       NAMPT   ITGB1    NAMPT_ITGB1_M_Monocyte_CD14_resting_L_T_Proliferating                       1.31                   5
##  6 M     L_T_TIM3._CD38._HLADR.  M_Monocyte_CD14_resting IFNG    IFNGR2   IFNG_IFNGR2_L_T_TIM3._CD38._HLADR._M_Monocyte_CD14_resting                  1.30                   6
##  7 M     L_T_TIM3._CD38._HLADR.  M_Monocyte_CD16         IFNG    IFNGR1   IFNG_IFNGR1_L_T_TIM3._CD38._HLADR._M_Monocyte_CD16                          1.30                   7
##  8 S     M_Monocyte_CD14_resting M_Monocyte_CD14_resting HLA.DMA CD4      HLA.DMA_CD4_M_Monocyte_CD14_resting_M_Monocyte_CD14_resting                 1.29                   8
##  9 M     M_Monocyte_CD16         L_T_Proliferating       HLA.DRA LAG3     HLA.DRA_LAG3_M_Monocyte_CD16_L_T_Proliferating                              1.28                   9
## 10 M     L_T_Proliferating       M_Monocyte_CD14_resting IFNG    IFNGR1   IFNG_IFNGR1_L_T_Proliferating_M_Monocyte_CD14_resting                       1.28                  10
```

``` r
ligand_oi = "IFNG"
receptor_oi = "IFNGR1"
group_oi = "M"
sender_oi = "L_T_TIM3._CD38._HLADR."
receiver_oi = "M_Monocyte_CD16"
```

``` r
p_violin = make_ligand_receptor_violin_plot(sce_sender = sce, sce_receiver = sce, ligand_oi = ligand_oi, receptor_oi = receptor_oi, group_oi = group_oi, group_id = group_id, sender_oi = sender_oi, receiver_oi = receiver_oi, sample_id = sample_id, celltype_id_sender = celltype_id, celltype_id_receiver = celltype_id)
p_violin
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

For `make_ligand_receptor_feature_plot`, your SingleCellExperiment
object should have a dimensionality reduction element stored.

## Zoom in on specific ligand-target interactions: show their expression in the single-cell data!

Make target gene violin and feature plots: `make_target_violin_plot` and
`make_target_feature_plot`

LILRB1: interesting gene in MIS-C CD16+ monocytes

``` r
target_oi = "LILRB1"

make_target_violin_plot(sce_receiver = sce, target_oi = target_oi, receiver_oi = receiver_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id_receiver = celltype_id)
```

![](detailed_analysis_steps_MISC_files/figure-gfm/unnamed-chunk-78-1.png)<!-- -->
