#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

#Get a list of projects
gdcprojects = getGDCprojects()
getProjectSummary('TCGA-BRCA')

luad_cnv_list <- c("TCGA-50-5932-01A-11D-1752-01", "TCGA-49-6742-01A-11D-1854-01",
                  "TCGA-44-6147-01A-11D-A273-01", "TCGA-55-6979-01A-11D-1943-01",
                  "TCGA-50-5931-01A-11D-1752-01", "TCGA-38-4626-01A-01D-1549-01",
                  "TCGA-91-6835-01A-11D-1854-01", "TCGA-49-4490-01A-21D-1854-01",
                  "TCGA-55-6970-01A-11D-1943-01", "TCGA-55-6969-01A-11D-1943-01",
                  "TCGA-91-6831-01A-11D-1854-01", "TCGA-50-5936-01A-11D-1624-01",
                  "TCGA-44-2665-01A-01D-0944-01", "TCGA-44-6776-01A-11D-1854-01",
                  "TCGA-55-6978-01A-11D-1943-01", "TCGA-49-6744-01A-11D-1854-01",
                  "TCGA-55-6982-01A-11D-1943-01", "TCGA-44-3396-01A-01D-1204-01",
                  "TCGA-49-6761-01A-31D-1943-01", "TCGA-44-6778-01A-11D-1854-01",
                  "TCGA-44-6145-01A-11D-1752-01", "TCGA-50-6595-01A-12D-1854-01",
                  "TCGA-44-6148-01A-11D-1752-01", "TCGA-91-6828-01A-11D-1854-01",
                  "TCGA-91-6849-01A-11D-1943-01", "TCGA-91-6847-01A-11D-1943-01",
                  "TCGA-44-2668-01A-01D-0944-01", "TCGA-44-2662-01A-01D-0944-01",
                  "TCGA-44-5645-01A-01D-1624-01", "TCGA-55-6985-01A-11D-1943-01",
                  "TCGA-44-2657-01A-01D-1549-01", "TCGA-73-4676-01A-01D-1752-01",
                  "TCGA-50-5939-01A-11D-1624-01", "TCGA-50-5935-01A-11D-1752-01",
                  "TCGA-55-6981-01A-11D-1943-01", "TCGA-91-6829-01A-21D-1854-01",
                  "TCGA-44-6777-01A-11D-1854-01", "TCGA-55-6983-01A-11D-1943-01",
                  "TCGA-91-6836-01A-21D-1854-01", "TCGA-55-6986-01A-11D-1943-01",
                  "TCGA-38-4625-01A-01D-1204-01", "TCGA-50-5933-01A-11D-1752-01",
                  "TCGA-38-4632-01A-01D-1752-01", "TCGA-55-6971-01A-11D-1943-01",
                  "TCGA-55-6975-01A-11D-1943-01", "TCGA-49-6743-01A-11D-1854-01",
                  "TCGA-49-4512-01A-21D-1854-01", "TCGA-55-6972-01A-11D-1943-01",
                  "TCGA-38-4627-01A-01D-1204-01", "TCGA-44-3398-01A-01D-1549-01",
                  "TCGA-44-2655-01A-01D-0944-01", "TCGA-44-6146-01A-11D-1752-01",
                  "TCGA-44-2661-01A-01D-1549-01", "TCGA-49-6745-01A-11D-1854-01")

luad_rna_tumor <- c("TCGA-50-5932-01A-11R-1755-07", "TCGA-49-6742-01A-11R-1858-07",
                    "TCGA-44-6147-01A-11R-1755-07", "TCGA-55-6979-01A-11R-1949-07",
                    "TCGA-50-5931-01A-11R-1755-07", "TCGA-38-4626-01A-01R-1206-07",
                    "TCGA-91-6835-01A-11R-1858-07", "TCGA-50-5930-01A-11R-1755-07",
                    "TCGA-49-4490-01A-21R-1858-07", "TCGA-55-6970-01A-11R-1949-07",
                    "TCGA-55-6969-01A-11R-1949-07", "TCGA-91-6831-01A-11R-1858-07",
                    "TCGA-50-5936-01A-11R-1628-07", "TCGA-44-2665-01A-01R-0946-07",
                    "TCGA-44-6776-01A-11R-1858-07", "TCGA-55-6978-01A-11R-1949-07",
                    "TCGA-49-6744-01A-11R-1858-07", "TCGA-55-6982-01A-11R-1949-07",
                    "TCGA-44-3396-01A-01R-1206-07", "TCGA-49-6761-01A-31R-1949-07",
                    "TCGA-44-6778-01A-11R-1858-07", "TCGA-44-6145-01A-11R-1755-07",
                    "TCGA-50-6595-01A-12R-1858-07", "TCGA-44-6148-01A-11R-1755-07",
                    "TCGA-91-6828-01A-11R-1858-07", "TCGA-91-6849-01A-11R-1949-07",
                    "TCGA-91-6847-01A-11R-1949-07", "TCGA-44-2668-01A-01R-A278-07",
                    "TCGA-44-2662-01A-01R-0946-07", "TCGA-44-5645-01A-01R-1628-07",
                    "TCGA-55-6985-01A-11R-1949-07", "TCGA-55-6984-01A-11R-1949-07",
                    "TCGA-44-2657-01A-01R-1107-07", "TCGA-73-4676-01A-01R-1755-07",
                    "TCGA-50-5939-01A-11R-1628-07", "TCGA-50-5935-01A-11R-1755-07",
                    "TCGA-55-6981-01A-11R-1949-07", "TCGA-91-6829-01A-21R-1858-07",
                    "TCGA-55-6980-01A-11R-1949-07", "TCGA-44-6777-01A-11R-1858-07",
                    "TCGA-55-6983-01A-11R-1949-07", "TCGA-91-6836-01A-21R-1858-07",
                    "TCGA-55-6986-01A-11R-1949-07", "TCGA-38-4625-01A-01R-1206-07",
                    "TCGA-50-5933-01A-11R-1755-07", "TCGA-38-4632-01A-01R-1755-07",
                    "TCGA-55-6971-01A-11R-1949-07", "TCGA-55-6975-01A-11R-1949-07",
                    "TCGA-49-6743-01A-11R-1858-07", "TCGA-49-4512-01A-21R-1858-07",
                    "TCGA-55-6972-01A-11R-1949-07", "TCGA-38-4627-01A-01R-1206-07",
                    "TCGA-44-3398-01A-01R-1107-07", "TCGA-44-2655-01A-01R-0946-07",
                    "TCGA-44-6146-01A-11R-A278-07", "TCGA-44-2661-01A-01R-1107-07",
                    "TCGA-49-6745-01A-11R-1858-07")

##build a query to retrieve data
query_TCGA_cnv <- GDCquery(project = 'TCGA-LUAD',
        data.category = 'Copy Number Variation',
        #sample.type = "Primary Tumor",
        data.type = "Gene Level Copy Number",
        workflow.type = 'ASCAT3',
        barcode = luad_cnv_list)

output_query_TCGA <- getResults(query_TCGA_cnv)


#build a query to retrieve gene expression data
query_TCGA_rna <- GDCquery(project = 'TCGA-LUAD',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor",
                       access = 'open',
                       barcode = luad_rna_tumor)


brca_rna_tumor <- c("TCGA-E9-A1RH-01A-21R-A169-07", "TCGA-BH-A1ET-01A-11R-A137-07",
                    "TCGA-BH-A0HK-01A-11R-A056-07", "TCGA-BH-A0H5-01A-21R-A115-07",
                    "TCGA-AC-A2FM-01A-11R-A19W-07", "TCGA-BH-A0DK-01A-21R-A056-07",
                    "TCGA-BH-A1FC-01A-11R-A13Q-07", "TCGA-E9-A1N9-01A-11R-A14D-07",
                    "TCGA-GI-A2C8-01A-11R-A16F-07", "TCGA-BH-A18L-01A-32R-A12D-07",
                    "TCGA-BH-A0BS-01A-11R-A12P-07", "TCGA-BH-A0BQ-01A-21R-A115-07",
                    "TCGA-BH-A209-01A-11R-A157-07", "TCGA-BH-A1EV-01A-11R-A137-07",
                    "TCGA-BH-A1FR-01A-11R-A13Q-07", "TCGA-BH-A0BC-01A-22R-A084-07",
                    "TCGA-BH-A0BA-01A-11R-A056-07", "TCGA-BH-A18J-01A-11R-A12D-07",
                    "TCGA-BH-A203-01A-12R-A169-07", "TCGA-BH-A18M-01A-11R-A12D-07",
                    "TCGA-BH-A1FJ-01A-11R-A13Q-07", "TCGA-BH-A18R-01A-11R-A12D-07",
                    "TCGA-E9-A1ND-01A-11R-A144-07", "TCGA-A7-A0DB-01A-11R-A00Z-07",
                    "TCGA-E2-A158-01A-11R-A12D-07", "TCGA-BH-A1FU-01A-11R-A14D-07",
                    "TCGA-E9-A1RF-01A-11R-A157-07", "TCGA-E9-A1N5-01A-11R-A14D-07",
                    "TCGA-E2-A1L7-01A-11R-A144-07", "TCGA-E2-A1IG-01A-11R-A144-07",
                    "TCGA-AC-A2FF-01A-11R-A17B-07", "TCGA-BH-A1FE-01A-11R-A13Q-07",
                    "TCGA-BH-A0DT-01A-21R-A12D-07", "TCGA-BH-A1EW-01A-11R-A137-07",
                    "TCGA-A7-A0CE-01A-11R-A00Z-07", "TCGA-BH-A1F2-01A-31R-A13Q-07",
                    "TCGA-A7-A13G-01A-11R-A13Q-07", "TCGA-BH-A18V-01A-11R-A12D-07",
                    "TCGA-AC-A2FB-01A-11R-A17B-07", "TCGA-A7-A13F-01A-11R-A12P-07",
                    "TCGA-BH-A0E0-01A-11R-A056-07", "TCGA-BH-A1F0-01A-11R-A137-07",
                    "TCGA-BH-A0B3-01A-11R-A056-07", "TCGA-BH-A1EU-01A-11R-A137-07",
                    "TCGA-BH-A0DD-01A-31R-A12P-07", "TCGA-AC-A23H-01A-11R-A157-07",
                    "TCGA-BH-A0AY-01A-21R-A00Z-07", "TCGA-E9-A1RB-01A-11R-A157-07",
                    "TCGA-E9-A1NA-01A-11R-A144-07", "TCGA-BH-A1FD-01A-11R-A13Q-07",
                    "TCGA-BH-A0BM-01A-11R-A056-07", "TCGA-BH-A1F8-01A-11R-A13Q-07",
                    "TCGA-E9-A1RI-01A-11R-A169-07", "TCGA-BH-A1F6-01A-11R-A13Q-07",
                    "TCGA-E9-A1N6-01A-11R-A144-07", "TCGA-A7-A0CH-01A-21R-A00Z-07",
                    "TCGA-BH-A0H7-01A-13R-A056-07", "TCGA-BH-A18U-01A-21R-A12D-07",
                    "TCGA-BH-A18U-01A-21R-A12D-07", "TCGA-E2-A1BC-01A-11R-A12P-07",
                    "TCGA-BH-A1FN-01A-11R-A13Q-07", "TCGA-BH-A1FM-01A-11R-A13Q-07",
                    "TCGA-E2-A15M-01A-11R-A12D-07", "TCGA-BH-A0DH-01A-11R-A084-07",
                    "TCGA-BH-A0AU-01A-11R-A12P-07", "TCGA-E9-A1RD-01A-11R-A157-07",
                    "TCGA-BH-A0C0-01A-21R-A056-07", "TCGA-BH-A0HA-01A-11R-A12P-07",
                    "TCGA-BH-A0BZ-01A-31R-A12P-07", "TCGA-BH-A0BW-01A-11R-A115-07",
                    "TCGA-BH-A1EN-01A-11R-A13Q-07", "TCGA-BH-A1FB-01A-11R-A13Q-07",
                    "TCGA-E9-A1NF-01A-11R-A14D-07", "TCGA-BH-A0H9-01A-11R-A056-07",
                    "TCGA-BH-A0DP-01A-21R-A056-07", "TCGA-A7-A0DC-01A-11R-A00Z-07",
                    "TCGA-E2-A1LB-01A-11R-A144-07", "TCGA-E2-A15I-01A-21R-A137-07",
                    "TCGA-BH-A0BV-01A-11R-A00Z-07", "TCGA-BH-A18N-01A-11R-A12D-07",
                    "TCGA-E9-A1N4-01A-11R-A14M-07", "TCGA-BH-A0DQ-01A-11R-A084-07",
                    "TCGA-BH-A0BT-01A-11R-A12P-07", "TCGA-BH-A0B7-01A-12R-A115-07",
                    "TCGA-BH-A18K-01A-11R-A12D-07", "TCGA-A7-A0D9-01A-31R-A056-07",
                    "TCGA-E9-A1NG-01A-21R-A14M-07", "TCGA-BH-A0DO-01B-11R-A12D-07",
                    "TCGA-BH-A0B8-01A-21R-A056-07", "TCGA-BH-A1FG-01A-11R-A13Q-07",
                    "TCGA-BH-A18Q-01A-12R-A12D-07", "TCGA-E9-A1R7-01A-11R-A14M-07",
                    "TCGA-BH-A0DG-01A-21R-A12P-07", "TCGA-BH-A208-01A-11R-A157-07",
                    "TCGA-BH-A0DL-01A-11R-A115-07", "TCGA-BH-A0C3-01A-21R-A12P-07",
                    "TCGA-A7-A13E-01A-11R-A277-07", "TCGA-BH-A0E1-01A-11R-A056-07",
                    "TCGA-BH-A0B5-01A-11R-A12P-07", "TCGA-BH-A204-01A-11R-A157-07",
                    "TCGA-BH-A0AZ-01A-21R-A12P-07", "TCGA-E2-A15K-01A-11R-A12P-07",
                    "TCGA-BH-A1FH-01A-12R-A13Q-07", "TCGA-BH-A18P-01A-11R-A12D-07",
                    "TCGA-BH-A0BJ-01A-11R-A056-07", "TCGA-E2-A1LH-01A-11R-A14D-07",
                    "TCGA-E2-A153-01A-12R-A12D-07", "TCGA-BH-A18S-01A-11R-A12D-07",
                    "TCGA-BH-A1EO-01A-11R-A137-07", "TCGA-E9-A1RC-01A-11R-A157-07",
                    "TCGA-BH-A0DZ-01A-11R-A00Z-07", "TCGA-E2-A1LS-01A-12R-A157-07",
                    "TCGA-BH-A0DV-01A-21R-A12P-07")

brca_cnv_list <- c("TCGA-E9-A1RH-01A-21D-A166-01", "TCGA-BH-A1ET-01A-11D-A134-01",
                   "TCGA-BH-A0HK-01A-11D-A059-01", "TCGA-BH-A0H5-01A-21D-A111-01",
                   "TCGA-AC-A2FM-01A-11D-A19X-01", "TCGA-BH-A0DK-01A-21D-A059-01",
                   "TCGA-BH-A1FC-01A-11D-A13J-01", "TCGA-E9-A1N9-01A-11D-A14F-01",
                   "TCGA-GI-A2C8-01A-11D-A16C-01", "TCGA-BH-A18L-01A-32D-A12A-01",
                   "TCGA-BH-A0BS-01A-11D-A12N-01", "TCGA-BH-A0BQ-01A-21D-A111-01",
                   "TCGA-BH-A209-01A-11D-A160-01", "TCGA-BH-A1EV-01A-11D-A134-01",
                   "TCGA-BH-A1FR-01A-11D-A13J-01", "TCGA-BH-A0BC-01A-22D-A087-01",
                   "TCGA-BH-A0BA-01A-11D-A059-01", "TCGA-BH-A18J-01A-11D-A12A-01",
                   "TCGA-BH-A203-01A-12D-A166-01", "TCGA-BH-A18M-01A-11D-A12A-01",
                   "TCGA-BH-A1FJ-01A-11D-A13J-01", "TCGA-BH-A18R-01A-11D-A12A-01",
                   "TCGA-E9-A1ND-01A-11D-A141-01", "TCGA-A7-A0DB-01A-11D-A011-01",
                   "TCGA-E2-A158-01A-11D-A12A-01", "TCGA-BH-A1FU-01A-11D-A14F-01",
                   "TCGA-E9-A1RF-01A-11D-A160-01", "TCGA-E9-A1N5-01A-11D-A14F-01",
                   "TCGA-E2-A1L7-01A-11D-A141-01", "TCGA-E2-A1IG-01A-11D-A141-01",
                   "TCGA-AC-A2FF-01A-11D-A17C-01", "TCGA-BH-A1FE-01A-11D-A13J-01",
                   "TCGA-BH-A0DT-01A-21D-A12A-01", "TCGA-A7-A0CE-01A-11D-A011-01",
                   "TCGA-BH-A1F2-01A-31D-A13J-01", "TCGA-A7-A13G-01A-11D-A13J-01",
                   "TCGA-BH-A18V-01A-11D-A12A-01", "TCGA-AC-A2FB-01A-11D-A17C-01",
                   "TCGA-A7-A13F-01A-11D-A12N-01", "TCGA-BH-A0E0-01A-11D-A059-01",
                   "TCGA-BH-A1F0-01A-11D-A134-01", "TCGA-BH-A0B3-01A-11D-A059-01",
                   "TCGA-BH-A1EU-01A-11D-A134-01", "TCGA-BH-A0DD-01A-31D-A12N-01",
                   "TCGA-AC-A23H-01A-11D-A160-01", "TCGA-BH-A0AY-01A-21D-A011-01",
                   "TCGA-E9-A1RB-01A-11D-A160-01", "TCGA-E9-A1NA-01A-11D-A141-01",
                   "TCGA-BH-A1FD-01A-11D-A13J-01", "TCGA-BH-A0BM-01A-11D-A059-01",
                   "TCGA-BH-A1F8-01A-11D-A13J-01", "TCGA-E9-A1RI-01A-11D-A166-01",
                   "TCGA-BH-A1F6-01A-11D-A13J-01", "TCGA-E9-A1N6-01A-11D-A141-01",
                   "TCGA-A7-A0CH-01A-21D-A011-01", "TCGA-BH-A0H7-01A-13D-A059-01",
                   "TCGA-BH-A18U-01A-21D-A12A-01", "TCGA-GI-A2C9-01A-11D-A21P-01",
                   "TCGA-E2-A1BC-01A-11D-A12N-01", "TCGA-BH-A1FN-01A-11D-A13J-01",
                   "TCGA-BH-A1FM-01A-11D-A13J-01", "TCGA-E2-A15M-01A-11D-A12A-01",
                   "TCGA-BH-A0DH-01A-11D-A087-01", "TCGA-BH-A0AU-01A-11D-A12N-01",
                   "TCGA-E9-A1RD-01A-11D-A160-01", "TCGA-BH-A0C0-01A-21D-A059-01",
                   "TCGA-BH-A0HA-01A-11D-A12N-01", "TCGA-BH-A0BZ-01A-31D-A12N-01",
                   "TCGA-BH-A0BW-01A-11D-A111-01", "TCGA-BH-A1EN-01A-11D-A13J-01",
                   "TCGA-BH-A1FB-01A-11D-A13J-01", "TCGA-E9-A1NF-01A-11D-A14F-01",
                   "TCGA-BH-A0H9-01A-11D-A059-01", "TCGA-BH-A0DP-01A-21D-A059-01",
                   "TCGA-A7-A0DC-01A-11D-A011-01", "TCGA-E2-A1LB-01A-11D-A141-01",
                   "TCGA-E2-A15I-01A-21D-A134-01", "TCGA-BH-A0BV-01A-11D-A011-01",
                   "TCGA-BH-A18N-01A-11D-A12A-01", "TCGA-E9-A1N4-01A-11D-A14J-01",
                   "TCGA-BH-A0DQ-01A-11D-A087-01", "TCGA-BH-A0BT-01A-11D-A12N-01",
                   "TCGA-BH-A0B7-01A-12D-A111-01", "TCGA-BH-A18K-01A-11D-A12A-01",
                   "TCGA-A7-A0D9-01A-31D-A059-01", "TCGA-E9-A1NG-01A-21D-A14J-01",
                   "TCGA-BH-A0DO-01B-11D-A12A-01", "TCGA-BH-A0B8-01A-21D-A059-01",
                   "TCGA-BH-A1FG-01A-11D-A13J-01", "TCGA-BH-A18Q-01A-12D-A12A-01",
                   "TCGA-E9-A1R7-01A-11D-A14J-01", "TCGA-BH-A0DG-01A-21D-A12N-01",
                   "TCGA-BH-A208-01A-11D-A160-01", "TCGA-BH-A0DL-01A-11D-A111-01",
                   "TCGA-BH-A0C3-01A-21D-A12N-01", "TCGA-A7-A13E-01A-11D-A12N-01",
                   "TCGA-BH-A0E1-01A-11D-A059-01", "TCGA-BH-A0B5-01A-11D-A12N-01",
                   "TCGA-BH-A204-01A-11D-A160-01", "TCGA-BH-A0AZ-01A-21D-A12N-01",
                   "TCGA-E2-A15K-01A-11D-A12N-01", "TCGA-BH-A1FH-01A-12D-A13J-01",
                   "TCGA-BH-A0BJ-01A-11D-A059-01", "TCGA-E2-A1LH-01A-11D-A14F-01",
                   "TCGA-E2-A153-01A-12D-A12A-01", "TCGA-BH-A18S-01A-11D-A12A-01",
                   "TCGA-BH-A1EO-01A-11D-A134-01", "TCGA-E9-A1RC-01A-11D-A160-01",
                   "TCGA-BH-A0DZ-01A-11D-A011-01", "TCGA-BH-A0DV-01A-21D-A12N-01")
  

query_TCGA_normal <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       access = 'open',
                       sample.type = "Solid Tissue Normal")
                       #barcode = brca_rna_tumor

query_TCGA_cnv <- GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Copy Number Variation',
                           #sample.type = "Primary Tumor",
                           data.type = "Gene Level Copy Number",
                           workflow.type = 'ASCAT3',
                           barcode = brca_cnv_list)

##build a query to retrieve data
getResults(query_TCGA_cnv)

#download
GDCdownload(query_TCGA_cnv)

#prepare data
brca_cnv <- GDCprepare(query_TCGA_cnv, summarizedExperiment = TRUE)
brca_cnv_tumor <- assay(brca_cnv, 'copy_number', rownames = TRUE)

brca_rna_tumor <- GDCprepare(query_TCGA_tumor, summarizedExperiment = TRUE)
brca_rna_tum <- assay(brca_rna_tumor, 'unstranded', rownames = TRUE)

gene_name <- as.data.frame(brca_cnv@rowRanges@elementMetadata@listData[["gene_name"]]) 
colnames(gene_name)[1] <- "GeneID"
brca_cnv_tumor <- as.data.frame(brca_cnv_tumor)
brca_cnv_tumor <- cbind(gene_name, brca_cnv_tumor)
brca_cnv_tumor <- brca_cnv_tumor[!duplicated(brca_cnv_tumor$GeneID), ] %>% remove_rownames %>% column_to_rownames(var="GeneID")

save(brca_cnv_tumor, file = '~/model_data/TCGA/breast_cancer/brca_cnv_tumor.Rdata')

#common_cols <- intersect(colnames(luad_rna_normal), colnames(luad_rna_tumor))
#luad_rna_normal <- luad_rna_normal[, c(-8,-14,-28,-30,-32,-33,-39,-40,-42,-47,-55,-57,-58)]

colnames(luad_rna_normal) <- c('s1_normal','s2_normal','s3_normal','s4_normal','s5_normal', 's6_normal',
                              's7_normal','s8_normal','s9_normal','s10_normal','s11_normal', 's12_normal',
                              's13_normal','s14_normal','s15_normal','s16_normal','s17_normal', 's18_normal',
                              's19_normal','s20_normal','s21_normal','s22_normal','s23_normal', 's24_normal',
                              's25_normal','s26_normal','s27_normal','s28_normal','s29_normal', 's30_normal',
                              's31_normal','s32_normal','s33_normal','s34_normal','s35_normal', 's36_normal',
                              's37_normal','s38_normal','s39_normal','s40_normal','s41_normal', 's42_normal',
                              's43_normal','s44_normal','s45_normal','s46_normal')
#rownames(countdata) <- seqdata[,1]
#substr("ThisIsAString", start=1, stop=5) #shorten sample names
#colnames(countdata) <- substr(colnames(countdata), 1, 7)
#colnames(countdata)==sampleinfo$SampleName
#all(colnames(countdata)==sampleinfo$SampleName)

#Filtering low counts genes
rna_normal_tumor_3 <- rna_normal_tumor_3[which(rowSums(rna_normal_tumor_3)>10),]
cnv <- cnv[(rownames(cnv) %in% rownames(rna_normal_tumor_3)),] #delete rows by name

save(luad_cnv_tumor, file = "~/model_data/TCGA/lung_cancer/cnv_luad.Rdata")
write.csv(rna_cnv, "model_data/TCGA/lung_cancer/last_test/rna_cnv.csv")

luad_cnv_tumor <- luad_cnv_tumor %>% remove_rownames %>% column_to_rownames(var="GeneID")
rna_normal <- rna_normal_tumor %>% select(1:46)
rna_tumor <- rna_normal_tumor %>% select(47:92)
cnv_tumor_3 <- luad_cnv_tumor %>% select(8,16,27)
rna_normal_tumor_3 <- rna_normal_tumor %>% as.data.frame() %>% select(8,16,27,54,62,73)
rna_normal_3 <- rna_normalized_normal %>% as.data.frame() %>% select(8,16,27)
rna_3 <- cbind(rna_tumor_3, rna_normal_3)

