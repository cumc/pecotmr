# requires ctwas package.
# remotes::install_github("xinhe-lab/ctwas",ref = "main")
#
#' Causal inference for TWAS with Summary Statistics
#'
#' @title cTWAS (causal-TWAS) wrapper
#' 
#' @description 
#' 
#' @param weight_matrices a list of weight vectors, one weight vector corresponds
#' to one gene, name of each list is gene name.
#' 
#' @param gwas_sumstat a data frame of gwas summary statistics with the columns of:
#' "id", "A1", "A2", "z". 
#' 
#' @param ld_dir a string, pointing to a directory containing all LD matrix files 
#' (*.cor.xz) and variant information (*.bim). The .bim file has 5 required columns: 
#' "chrom", "id", "pos", "alt", "ref" . 
#'
#' @param wd is a path as working directory that will later create folders for output
#' files and formatted LDs. 
#'
#' @param outputdir a string, the directory to store output
#' 
#' @param harmonize_z TRUE/FALSE. If TRUE, GWAS and eQTL genotype alleles are 
#' harmonized
#'
#' @param harmonize_wgt TRUE/FALSE. If TRUE, GWAS and eQTL genotype alleles are 
#' harmonized
#'
#' @param recover_strand_ambig_wgt TRUE/FALSE. If TRUE, a procedure is used to 
#' recover strand ambiguous variants. If FALSE, these variants are dropped from 
#' the prediction model. 

#' @param strand_ambig_action_z the action to take to harmonize strand ambiguous variants 
#' (A/T, G/C) between the z scores and LD reference. "drop" removes the ambiguous variant 
#' from the z scores. "none" treats the variant as unambiguous, flipping the z score to 
#' match the LD reference and then taking no additional action. "recover" imputes the sign 
#' of ambiguous z scores using unambiguous z scores and the LD reference and flips the z 
#' scores if there is a mismatch between the imputed sign and the observed sign of the z 
#' score. This option is computationally intensive.
#'
#' @param outname a string, the output name. 
#'
#' @return 
#'
#' 
#'
#'
# @importFrom ctwas impute_expr_z ctwas_rss
# @importFrom biomaRt useMart getBM
# @importFrom RSQLite dbConnect dbDisconnect dbWriteTable
#' @importFrom data.table fread read_delim
#' @importFrom readr read_delim
#' @importFrom reshape2 melt 
#' @importFrom purrr map

# work in progress 
# library("biomaRt")
# library("DBI")
# library("RSQLite")
# library("data.table")
# library("ctwas")
# library("reshape2")
# library("readr")
# library("map")
#
## raw example 
# z_snp <- data.table::fread("/mnt/vast/hpc/csg/cl4215/mrmash/workflow/ctwas_gwas.tsv", 
#     data.table=FALSE, header=TRUE)
#
# wgt_all <- readRDS("/mnt/vast/hpc/csg/cl4215/mrmash/workflow/weight_all.rds")
#
# region_list <- read.table("/mnt/vast/hpc/csg/snuc_pseudo_bulk/eight_celltypes_analysis/output/data_preprocessing/ALL/phenotype_data/ALL.log2cpm.region_list", 
#                          sep="\t", header=FALSE)
# colnames(region_list) <- c("chrom", "start", "stop", "id", "gene_name")
#
#
# k <- ctwas_wrapper(weight_matrices=wgt_all, 
#               gwas_sumstat=z_snp,
#               ld_dir="/mnt/vast/hpc/csg/cl4215/mrmash/workflow/LD_WGS/",
#               wd="/mnt/vast/hpc/csg/cl4215/mrmash/workflow/",
#               region_list=region_list,
#               plot=FALSE,
#               harmonize_z = T,
#               harmonize_wgt = T,
#               outname="test",
#               strand_ambig_action_z = "none",
#               recover_strand_ambig_wgt = T)
# 

ctwas_wrapper <- function(weight_matrices, 
                          gwas_sumstat,
                          ld_dir,
                          wd,
                          region_list,
                          plot=FALSE,
                          harmonize_z = T,
                          harmonize_wgt = T,
                          outname="test",
                          strand_ambig_action_z = "none",
                          recover_strand_ambig_wgt = T,
                          estimate_group_prior = T,
                          estimate_group_prior_var = T,
                          ncore = 16
                          ){
    
    
    # FORMAT WEIGHTS - > make db file
    all_snps <- c()
    for ( i in names(wgt_all)){
        all_snps <- c(all_snps, names(wgt_all[[i]]))
        }
    length(all_snps)

    wgt_table <- do.call(rbind, lapply(names(wgt_all), get_wgt_table, wgt_list=wgt_all)) 
    extra <- get_extra(wgt_list=wgt_all, region_list)
    
    extra <-  extra[extra$gene_type=="protein_coding",,drop=F]
    wgt_table <- wgt_table[wgt_table$gene %in% extra$gene,]
    
    mydb <- dbConnect(RSQLite::SQLite(), paste0(wd, "/weights.db"))
    dbWriteTable(mydb, "weights", wgt_table, overwrite=TRUE)
    dbWriteTable(mydb,"extra", extra, overwrite=TRUE)                            
    dbDisconnect(mydb)
    
    # FORMAT LD
    ifelse(!dir.exists(file.path(paste0(wd), "LD_format")), dir.create(file.path(paste0(wd), "LD_format")), FALSE)
    ld_R_dir <- paste0(wd, "/LD_format/") # output of formatted ld
    
    # get original LD
    lds <- list.files(path=ld_dir, pattern="*.cor.xz.bim", all.files=TRUE, full.names=FALSE)
    lds <- strsplit(lds, "*.cor.xz.bim")

    map(lds, write_ld_matrix, ld_dir, ld_R_dir, all_snps)
    
    # get pair-wise covariance /LD for harmonization
    ls <- list.files(path=ld_R_dir, pattern="*.RDS", all.files=TRUE, full.names=FALSE)
    ls <- strsplit(ls, ".RDS")
    
    covariance <- do.call(rbind, lapply(ls ,function(x){
            myld <- readRDS(paste0(ld_R_dir, "/", x, ".RDS"))
            Rvar <- data.table::fread(paste0(ld_R_dir, x, ".Rvar"), header = T)
    
            rownames(myld) <- paste0("chr", Rvar$chrom, "_", sub("^.*?\\:", "", Rvar$id), "_b38")
            colnames(myld) <- rownames(myld)
    
            cov <- reshape2::melt(myld)
            cov <- data.frame(t(apply(cov, 1, sort)))
            cov <- cov[duplicated(cov[, 1 : 2], MARGIN = 1), ]
            cov <- cov[, 3 : 1]
            colnames(cov) <- c("RSID1", "RSID2", "VALUE")
            return(cov)}
       )
    )
    
    data.table::fwrite(covariance, paste0(wd, "/weights.txt.gz"), sep="\t", 
                   quote = FALSE, row.names=FALSE)
    data.table::fwrite(region_list, paste0(wd, "/region_list.txt"), sep="\t", 
                   quote = FALSE, row.names=FALSE)
    
    ifelse(!dir.exists(file.path(paste0(wd), "output")), dir.create(file.path(paste0(wd), "output")), FALSE)
    
    res <- impute_expr_z(z_snp = z_snp,
                     weight = paste0(wd, "/weights.db"), 
                     ld_R_dir = ld_R_dir, 
                     outputdir = paste0(wd, "/output/"), 
                     outname = outname,
                     harmonize_z = harmonize_z,
                     harmonize_wgt = harmonize_wgt,
                     strand_ambig_action_z = strand_ambig_action_z,
                     recover_strand_ambig_wgt = T)

     
    # pars <- ctwas_rss(z_gene=res$z_gene,
    #              z_snp = res$z_snp, 
    #              ld_R_dir = ld_R_dir, 
    #              ld_regions_custom = paste0(wd, "/region_list.txt"), 
    #              outputdir = paste0(wd, "/output/"), 
    #              outname = outname,
    #              ncore = ncore)

    return(res)
    
}



get_wgt_table <- function(gene, wgt_list){
    wgt_table <- data.frame(gene=gene, 
                     rsid=names(wgt_list[[gene]]), 
                     varID=paste0(gsub( ":.*$", "", names(wgt_list[[gene]])) , "_" ,sub("^.*?\\:", "", names(wgt_list[[gene]])),"_b38"),
                     ref_allele=gsub( "_.*$", "", sub("^.*?\\_", "", names(wgt_list[[gene]]))), 
                     eff_allele=sub("^.*?\\_", "", sub("^.*?\\_", "", names(wgt_list[[gene]]))), 
                     weight=wgt_list[[gene]])
    return(wgt_table)
}


get_extra <- function(wgt_list, region_list){
    mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    gtype_table <- getBM(values = names(wgt_all), attributes = c("ensembl_gene_id", "gene_biotype"), mart=mart)
    
    extra <- data.frame(
        gene=names(wgt_list),
        genename=unlist(lapply(names(wgt_list), function(x) region_list$gene_name[region_list$id == x])),
        n.snps.in.model=unlist(lapply(names(wgt_list), function(x) length(wgt_list[[x]]))),
        gene_type=gtype_table$gene_biotype[gtype_table$ensembl_gene_id %in% names(wgt_list)], 
        pred.perf.R2=NA,
        pred.perf.pval=NA,
        pred.perf.qval=NA
        )
    return(extra)  
}


# function to format LD matrices for cTWAS
write_ld_matrix <- function(ld_dir, ld_R_dir, all_snps) {
  rvar_file <- data.table::fread(paste0(ld_dir, "/", i, ".cor.xz.bim"), header = FALSE)
  colnames(rvar_file) <- c("chrom", "id", "posg", "pos", "alt", "ref")
  rvar_file <- rvar_file[, c("chrom", "id", "pos", "alt", "ref")]
  indx <- which(rvar_file$id %in% all_snps)

  if (length(indx) >= 2) {
    rvar <- rvar_file[indx, ]
    rvar$variance <- 1
    rvar$chrom <- as.integer(rvar$chrom)

    ld_fname <- paste0(ld_dir, "/", i, ".cor")
    if (!file.exists(ld_fname)) {
      system(paste("xz -dk", paste0(ld_fname, ".xz")))
    }
    ld <- data.table::fread(paste0(ld_dir, "/", i, ".cor"), sep = " ")
    ld <- as.matrix(ld)
    ld <- ld[indx, indx]
    ld[upper.tri(ld)] <- t(ld)[upper.tri(ld)]
    start <- min(rvar$pos)
    stop <- max(rvar$pos)
    CHR <- as.integer(unique(rvar$chrom))
    if (length(CHR) >= 2) {
      stop("has more than one chromosome provided in the LD file")
    }

    saveRDS(ld, paste0(ld_R_dir, "/ld_chr", CHR, ".R_snp.", start, "_", stop, ".RDS"))
    data.table::fwrite(rvar,
      file = paste0(ld_R_dir, "/ld_chr", CHR, ".R_snp.", start, "_", stop, ".Rvar"),
      sep = "\t", quote = F
    )
  }
}



# select variants for ctwas weights input
select_ctwas_weights <- function(weight_db_files, conditions = NULL,
                                 variable_name_obj = c("preset_variants_result", "variant_names"),
                                 susie_obj = c("preset_variants_result", "susie_result_trimmed"),
                                 twas_weights_table = "twas_weights", max_var_selection,
                                 min_rsq_threshold = 0.01, p_val_cutoff = 0.05) {
  ## Internal function to load and validate data from RDS files
  load_and_validate_data <- function(weight_db_files, conditions, variable_name_obj) {
    all_data <- lapply(weight_db_files, readRDS)
    unique_regions <- unique(unlist(lapply(all_data, function(data) names(data))))

    # Check if region from all RDS files are the same
    if (length(unique_regions) != 1) {
      stop("The RDS files do not refer to the same region.")
    } else {
      # Assuming all data refer to the same region, now combine data by conditions
      combined_all_data <- do.call("c", lapply(all_data, function(data) data[[1]]))
    }
    # Set default for 'conditions' if they are not specified
    if (is.null(conditions)) {
      conditions <- names(combined_all_data) # get only condition name
    }
    ## Check if the specified condition and variable_name_obj are available in all files
    if (!all(conditions %in% names(combined_all_data))) {
      stop("The specified condition is not available in all RDS files.")
    }
    return(combined_all_data)
  }
  # determine if the region is imputable and select the best model
  # Function to pick the best model based on adj_rsq and p-value
  pick_best_model <- function(combined_all_data, conditions, min_rsq_threshold, p_val_cutoff) {
    best_adj_rsq <- min_rsq_threshold
    # Extract performance table
    performance_tables <- lapply(conditions, function(condition) {
      get_nested_element(combined_all_data, c(condition, "twas_cv_result", "performance"))
    })
    names(performance_tables) <- conditions
    # Determine if a gene/region is imputable and select the best model
    model_selection <- lapply(conditions, function(condition) {
      best_model <- NULL
      available_models <- do.call(c, lapply(names(performance_tables[[condition]]), function(model) {
        if (!is.na(performance_tables[[condition]][[model]][, "adj_rsq"])) {
          return(model)
        }
      }))
      if (length(available_models) <= 0) {
        print(paste0("No model provided TWAS cross validation performance metrics information at condition ", condition, ". "))
        return(NULL)
      }
      for (model in available_models) {
        model_data <- performance_tables[[condition]][[model]]
        if (model_data[, "adj_rsq"] >= best_adj_rsq) {
          best_adj_rsq <- model_data[, "adj_rsq"]
          best_model <- model
        }
      }
      region_names <- do.call(c, lapply(conditions, function(condition) {
        combined_all_data[[condition]]$region_info$region_name
      }))
      if (is.null(best_model)) {
        print(paste0(
          "No model has p-value < ", p_val_cutoff, " and r2 >= ", min_rsq_threshold, ", skipping condition ", condition,
          " at genes/regions ", unique(region_names), ". "
        ))
        return(NULL) # No significant model found
      } else {
        best_model <- unlist(strsplit(best_model, "_performance"))
        print(paste0("The best model for condition ", condition, " at region ", unique(region_names), " is ", best_model, ". "))
        return(best_model)
      }
    })
    names(model_selection) <- conditions
    return(model_selection)
  }
  # Based on selected best model and imputable conditions, select variants based on susie output
  ctwas_select <- function(combined_all_data, conditions, imputable_status, max_var_selection, variable_name_obj, susie_obj) {
    names(imputable_status) <- conditions
    var_selection_list <- lapply(conditions, function(condition) {
      if (imputable_status[condition] == "non_imputable") {
        out <- list(
          variant_selection = rep(NA, max_var_selection),
          susie_result_trimmed = list(
            sets = list(cs = c(NULL)),
            pip = rep(NA, max_var_selection)
          )
        )
        return(out)
      }
      # if the condition is imputable
      out <- list(
        variant_names = get_nested_element(combined_all_data, c(condition, variable_name_obj)),
        susie_result_trimmed = get_nested_element(combined_all_data, c(condition, susie_obj))
      )
      # Go through cs to select top variants from each set until we select all variants
      select_cs_var <- function(set_cs_list, max_var_selection) {
        selected_indices <- c()
        while (length(selected_indices) < max_var_selection) {
          for (cs in set_cs_list$sets$cs) {
            if (length(selected_indices) < max_var_selection) {
              # Filter out any indices that have already been selected
              available_indices <- setdiff(cs, selected_indices)
              if (length(available_indices) > 0) {
                # Get the index with the highest PIP score
                top_index <- available_indices[which.max(set_cs_list$pip[available_indices])]
                selected_indices <- c(selected_indices, top_index)
              }
            }
          }
        }
        return(selected_indices)
      }
      # check if the condition has top_loci table
      if ("top_loci" %in% names(combined_all_data[[condition]][["preset_variants_result"]])) {
        out$top_loci <- get_nested_element(combined_all_data, c(condition, "preset_variants_result", "top_loci"))
        if (length(out$top_loci[, "variant_id"]) <= max_var_selection) {
          out$variant_selection <- out$top_loci$variant_id
        } else {
          # if top_loci variants more than max_var_selection, we loop through CS set to select top pip variants from each CS set
          # out$out$susie_result_trimmed is from "preset_variants_result"-"susie_result_trimmed"
          if (!is.null(out$susie_result_trimmed$sets$cs)) {
            # check total variant number in the CS sets
            totl_cs_indice <- do.call(sum, lapply(names(out$susie_result_trimmed$sets$cs), function(L) {
              length(out$susie_result_trimmed$sets$cs[[L]])
            }))
            if (max_var_selection <= totl_cs_indice) {
              selec_idx <- select_cs_var(out$susie_result_trimmed, max_var_selection)
              out$variant_selection <- out$variant_names[selec_idx]
            } else {
              # when top loci has larger number of variants than max_var_selection, but CS sets has less number than max_var_selection
              selec_idx <- unlist(out$susie_result_trimmed$sets$cs)
              out$variant_selection <- out$variant_names[selec_idx]
            }
          } else {
            top_idx <- order(out$susie_result_trimmed$pip, decreasing = TRUE)[1:max_var_selection]
            out$variant_selection <- out$variant_names[top_idx]
          }
        }
      } else {
        # if the condition did not come with the top_loci table, we select top pip variants from [[condition]]$preset_variants_result$susie_result_trimmed$pip
        print(paste0(condition, " do not have top_loci table, select top pip variants. "))
        top_idx <- order(out$susie_result_trimmed$pip, decreasing = TRUE)[1:max_var_selection]
        out$variant_selection <- out$variant_names[top_idx]
      }
      return(out)
    })
    names(var_selection_list) <- conditions
    return(var_selection_list)
  }

  # Internal function to align and merge weight matrices
  align_and_merge <- function(weights_list, variable_objs) {
    # Get the complete list of variant names across all files
    all_variants <- unique(unlist(variable_objs))
    consolidated_list <- list()
    # Fill the matrix with weights, aligning by variant names
    for (i in seq_along(weights_list)) {
      # Initialize the temp matrix with zeros
      existing_colnames <- character(0)
      temp_matrix <- matrix(0, nrow = length(all_variants), ncol = ncol(weights_list[[i]]))
      rownames(temp_matrix) <- all_variants
      idx <- match(variable_objs[[i]], all_variants)
      temp_matrix[idx, ] <- weights_list[[i]]
      # Ensure no duplicate column names
      new_colnames <- colnames(weights_list[[i]])
      dups <- duplicated(c(existing_colnames, new_colnames))
      if (any(dups)) {
        duplicated_names <- paste(c(existing_colnames, new_colnames)[dups], collapse = ", ")
        stop("Duplicate column names detected during merging process: ", duplicated_names, ".")
      }
      existing_colnames <- c(existing_colnames, new_colnames)
      consolidated_list[[i]] <- temp_matrix
      colnames(consolidated_list[[i]]) <- existing_colnames
    }
    return(consolidated_list)
  }

  # Internal function to consolidate weights for given condition
  consolidate_weights_list <- function(combined_all_data, conditions, variable_name_obj, twas_weights_table) {
    # Set default for 'conditions' if they are not specified
    if (is.null(conditions)) {
      conditions <- names(combined_all_data)
    }
    combined_weights_by_condition <- lapply(conditions, function(condition) {
      temp_list <- get_nested_element(combined_all_data, c(condition, twas_weights_table))
      temp_list <- temp_list[!names(temp_list) %in% "variant_names"]
      sapply(temp_list, cbind)
    })
    names(combined_weights_by_condition) <- conditions
    if (is.null(variable_name_obj)) {
      # Standard processing: Check for identical row numbers and consolidate
      row_numbers <- sapply(combined_weights_by_condition, function(data) nrow(data))
      if (length(unique(row_numbers)) > 1) {
        stop("Not all files have the same number of rows for the specified condition.")
      }
      weights <- combined_weights_by_condition
    } else {
      # Processing with variable_name_obj: Align and merge data, fill missing with zeros
      variable_objs <- lapply(conditions, function(condition) {
        get_nested_element(combined_all_data, c(condition, variable_name_obj))
      })
      weights <- align_and_merge(combined_weights_by_condition, variable_objs)
    }
    names(weights) <- conditions
    return(weights)
  }

  ## Load, validate, and consolidate data
  combined_all_data <- load_and_validate_data(weight_db_files, conditions, variable_name_obj)
  ## we first select best model, then select variants based on susie obj output, then extract weight value for these variants in selected model
  model_selection <- pick_best_model(combined_all_data, conditions, min_rsq_threshold, p_val_cutoff)
  region_info <- combined_all_data[[1]]$region_info
  if (is.null(unlist(model_selection))) {
    print("No model meets the p_value threshold and R-squared minimum in all conditions. Region is not imputable. ")
    model_selection$imputable <- FALSE
    return(list(susie_results = NULL, model_selection = model_selection, weights = NULL, region_info = region_info))
  } else {
    model_selection$imputable <- TRUE
  }
  # check for non-imputable conditions for this gene
  names(conditions) <- rep("impputable", length(conditions))
  non_imp_conditions <- na.omit(unlist(lapply(conditions, function(condition) {
    if (is.null(model_selection[[condition]])) {
      return(condition)
    }
  })))
  if (length(non_imp_conditions) >= 1) {
    names(conditions)[which(conditions %in% non_imp_conditions)] <- "non_imputable"
  }
  # at imputable conditions, sometimes susie_obj is not available.
  # if susie_obj is not available
  if (model_selection$imputable == TRUE) {
    for (k in conditions) {
      if (is.null(combined_all_data[[k]][[susie_obj[1]]][[susie_obj[2]]])) {
        model_selection[k] <- list(NULL)
        names(conditions)[match(k, conditions)] <- "non_imputable"
      }
    }
  }
  # run ctwas_select on all conditions first to keep list structure the same, regardless of imputable status.
  ctwas_select_result <- ctwas_select(combined_all_data, conditions, names(conditions), max_var_selection, variable_name_obj, susie_obj)
  weights <- consolidate_weights_list(combined_all_data, conditions, variable_name_obj, twas_weights_table)
  # assign the weights_selection as NULL in ctwas_select_result #
  if (length(non_imp_conditions) >= 1) {
    for (unimp_conds in non_imp_conditions) {
      ctwas_select_result[[unimp_conds]][["variant_selection"]] <- rep(NA, length(ctwas_select_result[[unimp_conds]][["variant_selection"]]))
    }
  }
  # return results
  return(list(
    susie_results = ctwas_select_result, model_selection = model_selection,
    weights = weights, region_info = region_info
  ))
}


# ctwas input formatting functions for multigroup ctwas
get_ctwas_input <- function(summary_report, xqtl_meta_data, outdir, outname, chr_input_file = NULL) {
  if (!is.null(chr_input_file)) {
    out <- readRDS(chr_input_file)
  } else {
    out <- NULL
  }
  # only process new genes and add to pre-existing ctwas input file
  summary_report <- summary_report[summary_report$IsImputable, ]
  genes <- summary_report$gene[summary_report$gene %in% xqtl_meta_data$region_id]
  if (!is.null(chr_input_file)) {
    processed_genes <- unique(do.call(c, lapply(names(out[[1]]), function(chr) {
      gsub("\\|.*$", "", names(out[[1]][[chr]]$weights))
    })))
    genes <- genes[!genes %in% processed_genes]
  }
  file_list <- na.omit(summary_report$path[summary_report$gene %in% genes])
  if (length(file_list) >= 1) {
    # loop through genes in a chromosome
    foreach(file = file_list) %do% {
      twas_rs <- readRDS(file) # file is the ctwas variant selection output
      gene <- summary_report$gene[summary_report$path == file]
      chr <- as.integer(summary_report$chrom[summary_report$path == file])
      data_type <- summary_report$type[summary_report$path == file]
      studies <- na.omit(names(twas_rs))

      for (study in studies) {
        # update colnames of gwas sumstat z_snp
        gwas <- twas_rs[[study]]$gwas_qced
        id_idx <- match("variant_id", colnames(gwas$sumstats))
        colnames(gwas$sumstats)[id_idx] <- "id"
        # merge gwas sumstats and weights info
        out[[study]][[paste0("chr", chr)]][["z_snp"]] <- rbind(out[[study]][[paste0("chr", chr)]][["z_snp"]], gwas$sumstats[, c("id", "A1", "A2", "z")])
        out[[study]][[paste0("chr", chr)]][["z_snp"]] <- out[[study]][[paste0("chr", chr)]][["z_snp"]][!duplicated(out[[study]][[paste0("chr", chr)]][["z_snp"]][, "id"]), ]
        out[[study]][[paste0("chr", chr)]][["weights_list"]] <- c(out[[study]][[paste0("chr", chr)]][["weights_list"]], format_ctwas_weights(twas_rs, study, data_type))
      }
    }
    if (!is.null(chr_input_file)) {
      if (is.null(outdir)) outdir <- dirname(chr_input_file)
      if (is.null(outname)) outname <- ""
      print(paste0("saving formated output as ", outdir, "/", strsplit(basename(chr_input_file), ".rds")[[1]], "_", outname, "_updated.rds. "))
      saveRDS(out, paste0(outdir, "/", strsplit(basename(chr_input_file), ".rds")[[1]], "_", outname, "_updated.rds"))
      return(out)
    } else {
      if (is.null(outdir)) outdir <- getwd()
      print(paste0("saving formated output as ", outdir, "/ctwas_input_", outname, ".rds. "))
      saveRDS(out, paste0(outdir, "/ctwas_input_", outname, ".rds"))
      return(out)
    }
  }
  print(paste0("No ctwas input file found. "))
  return(NULL)
}


# extract selected weights from a gene
format_ctwas_weights <- function(twas_rs, study, type) {
  # Register parallel backend to use multiple cores
  region_info <- twas_rs[[study]][["region_info"]]
  chrom <- as.integer(region_info$grange$chrom)
  gene <- unique(region_info$region_name)
  # check if this gene is imputable at study
  if (twas_rs[[study]]$model_selection[[1]]$imputable) {
    conditions <- names(twas_rs[[study]]$weights)
    imputable_conditions <- unlist(lapply(conditions, function(condition) {
      ifelse(is.null(twas_rs[[study]][["model_selection"]][[condition]][["method"]]), return(NULL), return(condition))
    }))
    weights_list <- lapply(imputable_conditions, function(condition) {
      # remove NA weights
      wgt <- twas_rs[[study]][["weights"]][[condition]][["selected_weights"]]
      wgt_idx <- !is.na(wgt)
      wgt <- wgt[wgt_idx]
      if (length(wgt) == 0) {
        return(NULL)
      }
      variants <- twas_rs[[study]][["weights"]][[condition]][["variant_selection"]][wgt_idx]
      var_idx <- na.omit(match(variants, paste0("chr", colnames(twas_rs[[study]]$gwas_qced$LD))))
      # return NULL if no variants found in common of weights and LD
      if (length(var_idx) == 0) {
        return(NULL)
      }
      variants <- paste0("chr", colnames(twas_rs[[study]]$gwas_qced$LD))[var_idx] # get final variants name
      # get variant's variance to scale weights
      variance_idx <- match(variants, paste0("chr", twas_rs[[study]]$gwas_qced$sumstats$variant_id))
      variance_snp <- twas_rs[[study]]$gwas_qced$variance[variance_idx]
      # scale weights
      wgt <- wgt * sqrt(variance_snp) # wgt <- wgt*sqrt(ld_snpinfo_wgt$variance[ld_snpinfo_wgt.idx])
      wgt <- as.matrix(wgt[variants])
      colnames(wgt) <- "weight"
      # get LD - correlation matrix for selected variance
      ld_wgt <- as.matrix(twas_rs[[study]]$gwas_qced$LD[var_idx, var_idx])
      colnames(ld_wgt) <- rownames(ld_wgt) <- variants
      # p0 and p1 represent the min max position of weights variants
      p0 <- min(sapply(variants, function(x) {
        as.integer(unlist(strsplit(x, ":"))[2])
      }))
      p1 <- max(sapply(variants, function(x) {
        as.integer(unlist(strsplit(x, ":"))[2])
      }))
      return(list(
        chrom = chrom,
        p0 = p0,
        p1 = p1,
        wgt = wgt,
        R_wgt = ld_wgt,
        gene_name = gene,
        weight_name = gsub("_.*$", "", condition),
        type = type,
        context = gsub("_.*$", "", condition),
        n_wgt = length(var_idx)
      ))
    })
    names(weights_list) <- paste0(gene, "|", gsub("_.*$", "", imputable_conditions))
    return(weights_list)
  } else {
    print(paste0("Gene ", gene, " is not imputable, skipping. "))
  }
}




# extract ld file name from ld meta file
ld_meta_to_list <- function(ld_meta_file) {
  ld_meta <- data.table::fread(ld_meta_file, data.table = FALSE, header = TRUE)
  ld_file_list <- gsub(",.*$", "", ld_meta$path)
  return(ld_file_list)
}

# convert to symmetrical full matrix, save into RDS/Rvar format
format_ctwas_ld <- function(ld_paths, outdir) {
  region_info <- data.frame()
  for (ld_path in ld_paths) {
    ld <- as.matrix(read.table(ld_path, sep = " "))
    ld[lower.tri(ld)] <- t(ld)[lower.tri(ld)]
    saveRDS(ld, paste0(outdir, "/LD_", basename(ld_path), ".RDS"))

    bim <- read.table(paste0(ld_path, ".bim"), header = FALSE)[, -c(3, 9)] # remove posg and number missing sample column
    colnames(bim) <- c("chrom", "id", "pos", "alt", "ref", "variance", "allele_freq")
    bim$id <- gsub("_", ":", bim$id)
    data.table::fwrite(bim, file = paste0(outdir, "/LD_", basename(ld_path), ".Rvar"), sep = "\t", quote = F)
    print(paste0("writing ld file ", outdir, "/LD_", basename(ld_path), ".RDS"))
    print(paste0("writing ld file ", outdir, "/LD_", basename(ld_path), ".Rvar"))
    ld_info <- data.frame(
      chrom = as.integer(sub(".*chr", "", gsub("_.*$", "", basename(ld_path)))),
      start = as.integer(unlist(strsplit(basename(ld_path), "_"))[2]),
      stop = as.integer(gsub("\\..*", "", gsub(".*\\_", "", basename(ld_path)))),
      region_tag = paste0(gsub("_.*$", "", basename(ld_path)), ":", unlist(strsplit(basename(ld_path), "_"))[2], "-", gsub("\\..*", "", gsub(".*\\_", "", basename(ld_path)))),
      LD_matrix = paste0(outdir, "/LD_", basename(ld_path), ".RDS"),
      SNP_info = paste0(outdir, "/LD_", basename(ld_path), ".Rvar")
    )
    region_info <- rbind(region_info, ld_info)
  }
  return(region_info)
}


# function to scan all the ld files in a directory and generate region_info table
get_dir_region_info <- function(ld_dir) {
  ld_file_names <- list.files(path = ld_dir, pattern = "*.RDS", all.files = FALSE, full.names = FALSE)

  # Initialize vectors to store the positions, a2, and a1 elements
  chrom <- start <- stop <- vector("numeric", length(ld_file_names))
  # Extract the elements
  for (i in seq_along(ld_file_names)) {
    parts <- strsplit(ld_file_names[i], "_")[[1]]
    chrom[i] <- parts[2]
    start[i] <- parts[3]
    stop[i] <- gsub("\\..*", "", parts[4])
  }
  chrom <- as.integer(readr::parse_number(chrom))
  start <- as.integer(start)
  stop <- as.integer(stop)
  return(data.frame(
    chrom = chrom, start = start, stop = stop,
    region_tag = paste0("chr", chrom, ":", start, "-", stop),
    LD_matrix = paste0(ld_dir, "/", ld_file_names),
    SNP_info = paste0(ld_dir, "/", sub(".[^.]+$", "", ld_file_names), ".Rvar")
  ))
}
