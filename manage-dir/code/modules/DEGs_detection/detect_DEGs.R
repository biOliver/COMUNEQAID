library(foreach)
library(doParallel)
library(lme4)

run_glm = function(exp, metadata, model, mixed_model = T, family = gaussian()){
  tmp_data = data.frame(exp = exp, metadata)
  if (mixed_model) {
    fit_glm = glmer(model, data = tmp_data, family = family) # Using a mixed effect model.
  } else {
    fit_glm = glm(model, data = tmp_data, family = family) # Using a fixed effect model (typical glm; Super fast!!)
  }
  fit_glm.coef = summary(fit_glm)$coefficients # Return coefficients for all variables.
  return(fit_glm.coef)
}

snRNA_DEGs = function(object, fixed_effects, random_effects = NA, mixed_model = T, family = gaussian(), assay = 'SCT', slot = 'data', deg_type = 'abs', coefficient = 't value', n_cores = 50){
  registerDoParallel(n_cores)
  object.exp = GetAssayData(object, slot = slot, assay = assay)
  if (mixed_model & grepl('|', random_effects, fixed = T) == F) 
    return(message('Wrong random effect term. Use "|" to indicate random-effect terms, e.g., random_effects = (1|hash_id).'))
  if (mixed_model) {
    deg_model = paste(fixed_effects, random_effects, sep = '+')
  } else {
    deg_model = fix_effects
  }
  if (deg_type == 'rank') { 
    # Use rank test to ignore the magnitude of difference (might be useful to detect lowly expressed genes)
    tmp_model = as.formula(paste0('rank(exp) ~ ',deg_model))
  } else { 
    tmp_model = as.formula(paste0('exp ~ ',deg_model ))
  }
  target_genes = rownames(object.exp)
  exp_diff = data.frame(foreach(x = target_genes, .combine = rbind) %dopar% {
    # If there is an error, it will return 0. 
    # Only retain one coefficient (usually t value is most informative).
    tryCatch(run_glm(object.exp[x,], object@meta.data, tmp_model, mixed_model = mixed_model, family = family)[,coefficient], error = function(e) {return(0)})
  })
  rownames(exp_diff) = target_genes
  return(exp_diff)
}
