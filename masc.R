.sconline.MASCfn = function(
        dataset, 
        cluster, 
        contrast, 
        random_effects = NULL, 
        fixed_effects = NULL,
        verbose = TRUE, 
        jackknife=F,
        statistical.test="Wald") {
    
    
    #Adapted from Fonseka et al. PMID: 30333237
    
    # Check inputs
    require(lme4)
    if (is.factor(dataset[[contrast]]) == FALSE & is.numeric(dataset[[contrast]]) == FALSE) {
        stop("Specified contrast term should be coded as a factor or numeric in the dataset")
    }
    
    match.arg(statistical.test,c("LRT","Wald"))
    
    
    # Convert cluster assignments to string
    cluster = as.character(cluster)
    # Prepend design matrix generated from cluster assignments
    designmat = model.matrix(~ cluster + 0, data.frame(cluster = cluster))
    dataset = cbind(designmat, dataset)
    # Create output list to hold results
    res = vector(mode = "list", length = length(unique(cluster)))
    names(res) = attributes(designmat)$dimnames[[2]]
    
    # Create model formulas
    if (!is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs = paste0(c(paste0(fixed_effects, collapse = " + "),
                                                    paste0("(1|", random_effects, ")", collapse = " + ")),
                                                collapse = " + ")
        if (verbose == TRUE & statistical.test=="LRT") {
            message(paste("Using null model:", "cluster ~", model_rhs))
        }
    } else if (!is.null(fixed_effects) && is.null(random_effects)) {
        model_rhs = paste0(fixed_effects, collapse = " + ")
        if (verbose == TRUE&statistical.test=="LRT") {
            message(paste("Using null model:", "cluster ~", model_rhs))
            # For now, do not allow models without mixed effects terms
            
        }
        stop("No random effects specified")
    } else if (is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs = paste0("(1|", random_effects, ")", collapse = " + ")
        if (verbose == TRUE&statistical.test=="LRT") {
            message(paste("Using null model:", "cluster ~", model_rhs))
        }
    } else {
        model_rhs = "1" # only includes intercept
        if (verbose == TRUE&statistical.test=="LRT") {
            message(paste("Using null model:", "cluster ~", model_rhs))
            
        }
        stop("No random or fixed effects specified")
    }
    
    # Initialize list to store model objects for each cluster
    cluster_models = vector(mode = "list", length = length(attributes(designmat)$dimnames[[2]]))
    names(cluster_models) = attributes(designmat)$dimnames[[2]]
    
    # Run nested mixed-effects models for each cluster
    for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
        test_cluster = attributes(designmat)$dimnames[[2]][i]
        if (verbose == TRUE) {
            message(paste("Creating logistic mixed models for", test_cluster))
        }
        null_fm = as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                                                     model_rhs), collapse = ""))
        full_fm = as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                                                     model_rhs), collapse = ""))

        if (verbose == TRUE) {
            message(paste("Null model:", null_fm))
            message(paste("Full model:", full_fm))
        }
        
        # Run null and full mixed-effects models
        full_model = lme4::glmer(formula = full_fm, data = dataset,
                                                            family = binomial, nAGQ = 1, verbose = 0,
                                                            control = glmerControl(optimizer = "bobyqa"))
        
        
        # calculate confidence intervals for contrast term beta
        if(is.factor(dataset[[contrast]])){
            contrast_lvl2 = paste0(contrast, levels(dataset[[contrast]])[2])
            
        } else {
            contrast_lvl2 = contrast
        }
        
        contrast_ci = confint.merMod(full_model, method = "Wald",
                                                                    parm = contrast_lvl2)
        
        
        if(statistical.test=="Wald"){
            pval=summary(full_model)
            
            pval=pval$coefficients[contrast_lvl2,4]
        } else {
            null_model = lme4::glmer(
                formula = null_fm, data = dataset,
                family = binomial, nAGQ = 1, verbose = 0,
                control = glmerControl(optimizer = "bobyqa"))
            model_lrt = anova(null_model, full_model)
            pval=model_lrt[["Pr(>Chisq)"]][2]
        }
        
        # Save model objects to list
        cluster_models[[i]]$confint = contrast_ci
        cluster_models[[i]]$pval = pval
        cluster_models[[i]]$full_model = full_model
        
        #jackknifing
        jk_pvalvec=c()
        jk_coefvec=c()

        # jk_stable is a flag to indicate whether jackknifing was successful for all batches 
        # (NOT related to stability of the estimate sign or p value)
        # NOTE: CURRENTLY CAN ONLY JACKKNIFE WITH 1 RANDOM EFFECT
        jk_stable=1
        if(jackknife){
            njks_total = length(unique(dataset[,random_effects]))
            njks = 1
            for(ibatch in unique(dataset[[random_effects]])){
                print(paste("Running jackknife for batch",njks,"out of",njks_total))
                njks=njks+1
                tmp_dataset=dataset[dataset[[random_effects]]!=ibatch,]
                
                jk_full_model = tryCatch({
                    lme4::glmer(formula = full_fm, data = tmp_dataset,
                    family = binomial, nAGQ = 1, verbose = 0,
                    control = glmerControl(optimizer = "bobyqa"))},error=function(e){return(F)
                })
                
                # code below only runs if model above runs successfully
                if(class(jk_full_model)!=class(T)){
                    jk_coefvec=c(jk_coefvec,fixef(jk_full_model)[[contrast_lvl2]])
                    if(statistical.test=="Wald"){
                        tmp_pval=summary(jk_full_model)
                        tmp_pval=tmp_pval$coefficients[contrast_lvl2,4]
                        jk_pvalvec=c(jk_pvalvec,tmp_pval)
                    } else {
                        jk_null_model = tryCatch({
                            lme4::glmer(formula = null_fm, data = tmp_dataset,
                            family = binomial, nAGQ = 1, verbose = 0,
                            control = glmerControl(optimizer = "bobyqa"))},error=function(e) {return(F)})
                        
                        if(class(jk_null_model)!=class(T)){
                            jk_model_lrt = anova(jk_null_model, jk_full_model)
                            # calculate confidence intervals for contrast term beta
                            jk_pvalvec=c(jk_pvalvec,jk_model_lrt[["Pr(>Chisq)"]][2])
                        } else {
                            jk_stable=0
                        }
                    }
                    
                } else {
                    jk_stable=0
                }
                
            }
        } else {
            jk_pvalvec=(-1)
            jk_coefvec=(-1)
        }
        
        # cluster_models[[i]]$jk_pval_median = median(jk_pvalvec)
        # cluster_models[[i]]$jk_pval_mean = mean(jk_pvalvec)
        # cluster_models[[i]]$jk_pval_max = max(jk_pvalvec)
        
        # cluster_models[[i]]$jk_coef_median = median(jk_coefvec)
        # cluster_models[[i]]$jk_coef_mean = mean(jk_coefvec)
        # cluster_models[[i]]$jk_stable = jk_stable
        # cluster_models[[i]]$jk_coef_max = jk_coefvec[which(abs(jk_coefvec)==max(jk_coefvec))[1]]
        # cluster_models[[i]]$jk_coef_min = jk_coefvec[which(abs(jk_coefvec)==min(jk_coefvec))[1]]
    }
    
    # Organize results into output dataframe
    output = data.frame(cluster = attributes(designmat)$dimnames[[2]],
                                             size = colSums(designmat))
    output$model.pvalue = sapply(cluster_models, function(x) x$pval)
    output[[paste(contrast_lvl2, "OR", sep = ".")]] = sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
    output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] = sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
    output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] = sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))
    # output[[paste(contrast_lvl2,"JK","Max", "OR", sep = ".")]] = sapply(cluster_models, function(x) {if(x$jk_coef_max==(-1)){-1} else{exp(x$jk_coef_max)}})
    # output[[paste(contrast_lvl2,"JK","Min", "OR", sep = ".")]] = sapply(cluster_models, function(x) {if(x$jk_coef_min==(-1)){-1} else{exp(x$jk_coef_min)}})
    # output[[paste(contrast_lvl2,"JK","Mean", "OR", sep = ".")]] = sapply(cluster_models, function(x) {if(x$jk_coef_mean==(-1)){-1} else {exp(x$jk_coef_mean)}})
    # output[[paste(contrast_lvl2,"JK","Median", "OR", sep = ".")]] = sapply(cluster_models, function(x) {if(x$jk_coef_median==(-1)){-1} else {exp(x$jk_coef_median)}})
    # output[[paste(contrast_lvl2,"JK","Max", "pvalue", sep = ".")]] = sapply(cluster_models, function(x) {if(x$jk_pval_max==(-1)){-1} else {x$jk_pval_max}})
    # output[[paste(contrast_lvl2,"JK","Mean", "pvalue", sep = ".")]] = sapply(cluster_models, function(x) {if(x$jk_pval_mean==(-1)){-1} else {x$jk_pval_mean}})
    # output[[paste(contrast_lvl2,"JK","Median", "pvalue", sep = ".")]] = sapply(cluster_models, function(x) {if(x$jk_pval_median==(-1)){-1} else {x$jk_pval_median}})
    # output[[paste(contrast_lvl2,"JK","Stable", sep = ".")]] = sapply(cluster_models, function(x) x$jk_stable)
    
    return(output)
}

masc_fn_continuous = function(
    dataset, # metadata dataframe
    cluster_colname,
    test_colname,
    fixed_effects, 
    random_effects
){
    
    cluster = as.character(dataset[[cluster_colname]])
    
    # Prepend design matrix generated from cluster assignments
    designmat = model.matrix(~ cluster + 0, data.frame(cluster = cluster))
    dataset = cbind(designmat, dataset)

    res = vector(mode = "list", length = length(unique(cluster)))
    names(res) = attributes(designmat)$dimnames[[2]]

    model_rhs = paste0(c(paste0(c(test_colname, fixed_effects), collapse = " + "),
                            paste0("(1|", random_effects, ")", collapse = " + ")),
                            collapse = " + ")
    model_rhs

    null_rhs = paste0(c(paste0(fixed_effects, collapse = " + "),
                            paste0("(1|", random_effects, ")", collapse = " + ")),
                            collapse = " + ") 
    null_rhs                 

    # Initialize list to store model objects for each cluster
    cluster_models = vector(mode = "list",
                            length = length(attributes(designmat)$dimnames[[2]]))
    names(cluster_models) = attributes(designmat)$dimnames[[2]]

    models=list()
    null_models = list()
    pvals = list()
    # Run nested mixed-effects models for each cluster
    for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
        test_cluster = attributes(designmat)$dimnames[[2]][i]
        
        print(paste("Creating logistic mixed models for", test_cluster))
        
        full_fm = as.formula(paste0(c(paste0(test_cluster, " ~ "),
                                    model_rhs), collapse = ""))

        
        print(paste("Full model:", full_fm))
        
        
        # Run null and full mixed-effects models
        full_model = lme4::glmer(formula = full_fm, data = dataset,
                                family = binomial, nAGQ = 1, verbose = 0,
                                control = glmerControl(optimizer = "bobyqa"))


        null_fm = as.formula(paste0(c(paste0(test_cluster, " ~ "),
                                    null_rhs), collapse = ""))
        
        print(paste("null model:", null_fm))
    

        # Run null and full mixed-effects models
        null_model = lme4::glmer(formula = null_fm, data = dataset,
                                family = binomial, nAGQ = 1, verbose = 0,
                                control = glmerControl(optimizer = "bobyqa"))

        cluster = sub("cluster", "", names(full_model@frame)[[1]])

        model_lrt = anova(null_model, full_model)
        pval=model_lrt[["Pr(>Chisq)"]][2]

        models[[cluster]] = full_model
        null_models[[cluster]] = null_model
        pvals[[cluster]] = pval
    }

    for (m in models){
        cluster = sub("cluster", "", names(m@frame)[[1]])   
        models[[cluster]] = m
    }
    models = models[names(models) != ""]
    names(models)

    model_coef_list = list()
    for (name in names(models)){
        cat("\n", name)
        model = models[[name]]
        

        model_summary=summary(model)
        coef = model_summary$coefficients[test_colname, 1]
        pval=model_summary$coefficients[test_colname, 4]
        ci = confint.merMod(model, method = "Wald", parm = test_colname)
        

        model_coef_list[[name]] = list(
            cluster_name=name, 
            or=exp(coef), 
            pval=pval, 
            ci_lower=exp(ci[[1]]), 
            ci_upper = exp(ci[[2]])
        )
    }

    model_coef_df = as.data.frame(do.call(rbind, model_coef_list))
    model_coef_df$or = as.numeric(model_coef_df$or)
    model_coef_df$ci_lower = as.numeric(model_coef_df$ci_lower)
    model_coef_df$ci_upper = as.numeric(model_coef_df$ci_upper)

    model_coef_df$log2_or = log2(model_coef_df$or)
    model_coef_df$log2_ci_lower = log2(model_coef_df$ci_lower)
    model_coef_df$log2_ci_upper = log2(model_coef_df$ci_upper)

    model_coef_df$cluster_name = factor(model_coef_df$cluster_name, levels = model_coef_df$cluster_name[order(model_coef_df$or)])
    model_coef_df$pval = as.numeric(model_coef_df$pval)
    model_coef_df$rank = rank(model_coef_df$or)

    return(model_coef_df)

}