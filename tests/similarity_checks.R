for (R_file in dir(getwd(), pattern = "*.R$")) {
    source(file.path(getwd(), R_file))
}

taxa <- read.csv("Blueberry_genus_table.tsv", skip = 1, sep="\t")
rownames(taxa) <- taxa$X.OTU.ID
taxa$X.OTU.ID <- NULL

metadata <- read.csv("Blueberry_metadata.tsv", sep="\t")
rownames(metadata) <- metadata$sampleid
metadata$sampleid <- NULL

x <- Maaslin2(list("input_data" = taxa, "input_metadata" = metadata, "output" = "output", min_abundance = 0.1, min_prevalence = 0.02, min_variance = 0.01, max_significance = 0.2, normalization = "CSS"))
y <- Maaslin2(list("input_data" = taxa, "input_metadata" = metadata, "output" = "output", fixed_effects = "comparison", min_abundance = 0.1, min_prevalence = 0.02, min_variance = 0.01, max_significance = 0.2, normalization = "CSS"))
y <- Maaslin2(list("input_data" = taxa, "input_metadata" = metadata, "output" = "output", formula = "comparison", min_abundance = 0.1, min_prevalence = 0.02, min_variance = 0.01, max_significance = 0.2, normalization = "CSS"))



y <- Maaslin2::Maaslin2("input_data" = taxa, "input_metadata" = metadata, "output" = "output", min_abundance = 0.1, min_prevalence = 0.02, min_variance = 0.01, max_significance = 0.2, normalization = "CSS")

sum(x$fitted == y$fitted, na.rm=T)
sum(x$fitted != y$fitted, na.rm=T)

packrat:::recursivePackageDependencies("Maaslin2", ignore = "", lib.loc = .libPaths()[1])
packrat:::recursivePackageDependencies("Maaslin2Lite", ignore = "", lib.loc = .libPaths()[1])


taxa <- read.csv("HMP2_taxonomy.tsv", sep="\t")
rownames(taxa) <- taxa$ID
taxa$ID <- NULL

metadata <- read.csv("HMP2_metadata.tsv", sep="\t")
rownames(metadata) <- metadata$ID
metadata$ID <- NULL

taxa <- taxa[rownames(taxa) %in% rownames(metadata)[!is.na(metadata$age)],]
metadata <- metadata[rownames(metadata) %in% rownames(metadata)[!is.na(metadata$age)],]

metadata$interaction <- metadata$age * (metadata$dysbiosis_binary == "Yes")
metadata$age_squared <- metadata$age^2

x <- Maaslin2(list("input_data" = taxa, "input_metadata" = metadata, "output" = "output", fixed_effects = c("age", "age_squared"), min_abundance = 0.1, min_prevalence = 0.02, min_variance = 0.01, max_significance = 0.2, normalization = "TSS", plot_scatter = F, plot_heatmap = F, save_scatter = F, save_models = F, standardize=F))
y <- Maaslin2(list("input_data" = taxa, "input_metadata" = metadata, "output" = "output", formula = "poly(age, 2)", min_abundance = 0.1, min_prevalence = 0.02, min_variance = 0.01, max_significance = 0.2, normalization = "TSS", plot_scatter = F, plot_heatmap = F, save_scatter = F, save_models = F, standardize=F))

y <- Maaslin2::Maaslin2("input_data" = taxa, "input_metadata" = metadata, "output" = "output", fixed_effects = c("age", "diagnosis"), random_effects = c("subject"), reference = ('diagnosis,`CD`'), min_abundance = 0.1, min_prevalence = 0.02, min_variance = 0.01, max_significance = 0.2, normalization = "TSS", plot_scatter = F, plot_heatmap = F, save_scatter = F, save_models = F)

sum(x$fitted == y$fitted, na.rm=T)
sum(x$fitted != y$fitted, na.rm=T)













