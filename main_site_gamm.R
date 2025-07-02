# 加载必要的包
library(mgcv)
library(dplyr)
library(ggplot2)
library(readxl)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(openxlsx)
library(parallel)  # 并行计算包
library(doParallel)  # 并行计算包
library(foreach)     # 并行循环包
library(broom)       # 模型结果整理
library(moments)     # 偏度峰度计算

# 设置图形主题
theme_set(theme_minimal())

# ================================
# 路径设置
# ================================
# 设置输出路径
OUTPUT_PATH <- "C:/Users/Administrator/Desktop/3.Gonogo和Stroop比-主站点后测/output"

# 创建输出目录
if (!dir.exists(OUTPUT_PATH)) {
  dir.create(OUTPUT_PATH, recursive = TRUE)
  cat("创建output目录:", OUTPUT_PATH, "\n")
}

# ================================
# 数据采样选项设置
# ================================
# 设置为TRUE使用1%数据进行快速测试，FALSE使用完整数据
USE_SAMPLE_DATA <- TRUE  # 修改这里控制是否使用1%数据
SAMPLE_PERCENTAGE <- 1  # 采样比例，0.01表示1%

# 检测并设置并行计算核心数
n_cores <- detectCores() - 5  # 保留一个核心给系统
cat("检测到", detectCores(), "个CPU核心，将使用", n_cores, "个核心进行并行计算\n")

# 注册并行后端
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# ================================
# 函数定义
# ================================

# 并行Bootstrap函数
parallel_bootstrap <- function(data, formula_str, pred_data, subject_var, age_range, n_bootstrap) {
  
  cat("    使用", n_cores, "个核心并行计算Bootstrap置信区间 (", n_bootstrap, "samples)...\n")
  
  # 并行Bootstrap计算
  deriv_bootstrap <- foreach(i = 1:n_bootstrap, .combine = rbind,
                             .packages = c("mgcv"), 
                             .errorhandling = "remove") %dopar% {
                               
                               # Bootstrap抽样
                               boot_indices <- sample(nrow(data), replace = TRUE)
                               boot_data <- data[boot_indices, ]
                               
                               tryCatch({
                                 boot_model <- gam(as.formula(formula_str), 
                                                   data = boot_data,
                                                   method = "REML")
                                 
                                 boot_pred <- predict(boot_model, newdata = pred_data,
                                                      exclude = paste0("s(", subject_var, ")"))
                                 
                                 # 计算bootstrap导数
                                 boot_deriv_6mo <- c(0, diff(boot_pred) / diff(age_range))
                                 boot_deriv_6mo[1] <- boot_deriv_6mo[2]
                                 boot_deriv_annual <- boot_deriv_6mo * 2
                                 boot_deriv_increment <- boot_deriv_annual * 0.5
                                 
                                 return(boot_deriv_increment)
                                 
                               }, error = function(e) {
                                 return(rep(NA, length(age_range)))
                               })
                             }
  
  return(deriv_bootstrap)
}

# 倒U型模式检测函数
detect_inverted_u_pattern <- function(results, age_range,
                                      min_peak_age = 6, max_peak_age = 18) {
  
  predictions <- results$predictions
  derivatives <- results$derivatives
  significant <- results$significant
  
  # 寻找峰值点（导数从正变负）
  peak_candidates <- c()
  
  for (i in 2:(length(derivatives) - 1)) {
    if (!is.na(derivatives[i-1]) && !is.na(derivatives[i+1])) {
      if (derivatives[i-1] > 0 && derivatives[i+1] < 0 &&
          age_range[i] >= min_peak_age && age_range[i] <= max_peak_age) {
        peak_candidates <- c(peak_candidates, i)
      }
    }
  }
  
  if (length(peak_candidates) == 0) {
    return(list(
      has_inverted_u = FALSE,
      peak_age = NA,
      peak_value = NA,
      increase_period = c(NA, NA),
      decrease_period = c(NA, NA),
      increase_duration = 0,
      decrease_duration = 0
    ))
  }
  
  # 选择最高的峰值
  peak_values <- predictions[peak_candidates]
  best_peak_idx <- peak_candidates[which.max(peak_values)]
  peak_age <- age_range[best_peak_idx]
  peak_value <- predictions[best_peak_idx]
  
  # 确定上升期和下降期
  increase_mask <- significant & (age_range <= peak_age) & (derivatives > 0)
  decrease_mask <- significant & (age_range >= peak_age) & (derivatives < 0)
  
  increase_ages <- age_range[increase_mask]
  decrease_ages <- age_range[decrease_mask]
  
  increase_period <- if(length(increase_ages) > 0) c(min(increase_ages), max(increase_ages)) else c(NA, NA)
  decrease_period <- if(length(decrease_ages) > 0) c(min(decrease_ages), max(decrease_ages)) else c(NA, NA)
  
  increase_duration <- if(all(!is.na(increase_period))) diff(increase_period) else 0
  decrease_duration <- if(all(!is.na(decrease_period))) diff(decrease_period) else 0
  
  return(list(
    has_inverted_u = TRUE,
    peak_age = peak_age,
    peak_value = peak_value,
    increase_period = increase_period,
    decrease_period = decrease_period,
    increase_duration = increase_duration,
    decrease_duration = decrease_duration
  ))
}

# 发展阶段变化率计算函数
calculate_development_stage_rates <- function(results, age_range,
                                              stage_boundaries = c(6, 8, 9, 12, 13, 15, 16, 18)) {
  
  derivatives <- results$derivatives
  significant <- results$significant
  
  # 定义发展阶段
  stage_names <- c('early_rapid_6_8', 'middle_sustained_9_12',
                   'late_plateau_13_15', 'final_stable_16_18')
  stage_ranges <- list(
    c(6, 8),   # 早期快速发展阶段
    c(9, 12),  # 中期持续发展阶段
    c(13, 15), # 后期发展阶段
    c(16, 18)  # 稳定阶段
  )
  
  stage_results <- list()
  
  for (i in 1:length(stage_names)) {
    stage_name <- stage_names[i]
    stage_range <- stage_ranges[[i]]
    start_age <- stage_range[1]
    end_age <- stage_range[2]
    
    # 找到该阶段的年龄范围索引
    stage_mask <- (age_range >= start_age) & (age_range < end_age)
    stage_sig_mask <- stage_mask & significant
    
    if (any(stage_sig_mask)) {
      stage_derivatives <- derivatives[stage_sig_mask]
      stage_ages <- age_range[stage_sig_mask]
      
      # 计算统计量
      mean_rate <- mean(stage_derivatives, na.rm = TRUE)
      mean_abs_rate <- mean(abs(stage_derivatives), na.rm = TRUE)
      duration <- if(length(stage_ages) > 1) max(stage_ages) - min(stage_ages) else 0
      total_obs_period <- max(age_range) - min(age_range)
      duration_pct <- if(total_obs_period > 0) (duration / total_obs_period * 100) else 0
      age_range_str <- if(length(stage_ages) > 0) 
        paste0(round(min(stage_ages), 1), "-", round(max(stage_ages), 1)) else "无显著发展"
      
      stage_results[[paste0(stage_name, "_mean_rate")]] <- mean_rate
      stage_results[[paste0(stage_name, "_abs_mean_rate")]] <- mean_abs_rate
      stage_results[[paste0(stage_name, "_duration")]] <- duration
      stage_results[[paste0(stage_name, "_duration_pct")]] <- duration_pct
      stage_results[[paste0(stage_name, "_age_range")]] <- age_range_str
    } else {
      stage_results[[paste0(stage_name, "_mean_rate")]] <- 0
      stage_results[[paste0(stage_name, "_abs_mean_rate")]] <- 0
      stage_results[[paste0(stage_name, "_duration")]] <- 0
      stage_results[[paste0(stage_name, "_duration_pct")]] <- 0
      stage_results[[paste0(stage_name, "_age_range")]] <- "无显著发展"
    }
  }
  
  return(stage_results)
}

# 发展完成度计算函数
calculate_development_completion <- function(results, age_range,
                                             cutoff_ages = c(12, 13, 14, 15, 16)) {
  
  significant <- results$significant
  
  if (!any(significant)) {
    base_result <- list(
      sig_start = NA,
      sig_end = NA,
      sig_duration = 0
    )
    for (cutoff in cutoff_ages) {
      base_result[[paste0("completion_", cutoff)]] <- 0
      base_result[[paste0("pre_", cutoff, "_duration")]] <- 0
    }
    return(base_result)
  }
  
  # 找到显著发展期
  sig_ages <- age_range[significant]
  sig_start <- min(sig_ages)
  sig_end <- max(sig_ages)
  sig_duration <- sig_end - sig_start
  
  result <- list(
    sig_start = sig_start,
    sig_end = sig_end,
    sig_duration = sig_duration
  )
  
  # 计算每个cutoff_age的完成度
  for (cutoff_age in cutoff_ages) {
    pre_cutoff_end <- min(cutoff_age, sig_end)
    pre_cutoff_duration <- max(0, pre_cutoff_end - sig_start)
    completion_pct <- if(sig_duration > 0) (pre_cutoff_duration / sig_duration * 100) else 0
    
    result[[paste0("completion_", cutoff_age)]] <- completion_pct
    result[[paste0("pre_", cutoff_age, "_duration")]] <- pre_cutoff_duration
  }
  
  return(result)
}

# 计算GAM纵向模型（修改后版本）
compute_gam_longitudinal <- function(data, outcome_var, age_var = "age", 
                                     subject_var = "Subject", time_var = "timepoint",
                                     age_range = NULL, n_bootstrap = 3000, 
                                     k_basis = 10) {
  
  cat("    Fitting GAM model for", outcome_var, "...\n")
  cat("    使用修改后的显著性标准: derivatives_6mo_increment >=", significance_threshold, "\n")
  
  if (is.null(age_range)) {
    age_range <- seq(min(data[[age_var]], na.rm = TRUE), 
                     max(data[[age_var]], na.rm = TRUE), 
                     by = 0.1)
  }
  
  tryCatch({
    # 构建标准GAM模型（使用随机截距处理被试间差异）
    formula_str <- paste0(outcome_var, " ~ s(", age_var, ", k=", k_basis, ") + ", 
                          time_var, " + s(", subject_var, ", bs='re')")
    
    cat("    GAM Formula:", formula_str, "\n")
    
    gam_model <- gam(as.formula(formula_str), 
                     data = data,
                     method = "REML")
    
    # 模型摘要用于提取平滑项显著性
    model_summary <- summary(gam_model)
    
    # 提取平滑项显著性检验结果
    smooth_table <- model_summary$s.table
    age_smooth_row <- grep(paste0("s\\(", age_var, "\\)"), rownames(smooth_table))
    
    if (length(age_smooth_row) > 0) {
      smooth_f_value <- smooth_table[age_smooth_row[1], "F"]
      smooth_p_value <- smooth_table[age_smooth_row[1], "p-value"]
      smooth_ref_df <- smooth_table[age_smooth_row[1], "Ref.df"]
      smooth_edf <- smooth_table[age_smooth_row[1], "edf"]
    } else {
      smooth_f_value <- NA
      smooth_p_value <- NA
      smooth_ref_df <- NA
      smooth_edf <- sum(model_summary$edf)
    }
    
    # 模型诊断
    residuals <- residuals(gam_model)
    residual_std <- sd(residuals, na.rm = TRUE)
    
    # 生成预测数据
    # 使用第一个时间点和虚拟被试进行预测
    pred_data <- data.frame(
      age = age_range,
      timepoint = data[[time_var]][1],  # 使用第一个时间点
      Subject = "dummy_subject"  # 虚拟被试
    )
    names(pred_data)[names(pred_data) == "age"] <- age_var
    names(pred_data)[names(pred_data) == "timepoint"] <- time_var
    names(pred_data)[names(pred_data) == "Subject"] <- subject_var
    
    # 预测
    predictions <- predict(gam_model, newdata = pred_data, 
                           exclude = paste0("s(", subject_var, ")"))  # 排除随机效应
    pred_se <- predict(gam_model, newdata = pred_data, se.fit = TRUE,
                       exclude = paste0("s(", subject_var, ")"))$se.fit
    
    # 95%置信区间
    pred_lower <- predictions - 1.96 * pred_se
    pred_upper <- predictions + 1.96 * pred_se
    
    # 计算导数 (6个月的变化率)
    derivatives_6mo <- c(0, diff(predictions) / diff(age_range))
    derivatives_6mo[1] <- derivatives_6mo[2]  # 填充第一个值
    
    # 年化处理和6个月增量计算
    derivatives_annual <- derivatives_6mo * 2
    derivatives_6mo_increment <- derivatives_annual * 0.5
    
    cat("    计算年化增量: 6个月变化率 × 2 × 0.5 =", outcome_var, "的6个月纯增量\n")
    
    # 并行Bootstrap置信区间计算
    start_time <- Sys.time()
    deriv_bootstrap <- parallel_bootstrap(data, formula_str, pred_data, 
                                          subject_var, age_range, n_bootstrap)
    end_time <- Sys.time()
    
    cat("    并行Bootstrap完成，耗时:", round(as.numeric(end_time - start_time, units = "secs"), 2), "秒\n")
    
    # 计算置信区间
    valid_bootstrap <- !is.na(rowSums(deriv_bootstrap))
    if (sum(valid_bootstrap) > 10) {
      deriv_lower <- apply(deriv_bootstrap[valid_bootstrap, ], 2, 
                           function(x) quantile(x, 0.025, na.rm = TRUE))
      deriv_upper <- apply(deriv_bootstrap[valid_bootstrap, ], 2, 
                           function(x) quantile(x, 0.975, na.rm = TRUE))
    } else {
      deriv_std <- sd(derivatives_6mo_increment, na.rm = TRUE)
      if (is.na(deriv_std)) deriv_std <- 0.1
      deriv_lower <- derivatives_6mo_increment - 2 * deriv_std
      deriv_upper <- derivatives_6mo_increment + 2 * deriv_std
    }
    
    # 修改后的显著性检验：使用绝对值阈值
    significant <- abs(derivatives_6mo_increment) >= significance_threshold
    
    cat("    使用新标准检测到", sum(significant, na.rm = TRUE), "个显著发展点\n")
    
    # 连续性检验：显著区域需要连续至少3个点
    significant_cleaned <- rep(FALSE, length(significant))
    rle_result <- rle(significant)
    
    current_pos <- 1
    for (i in 1:length(rle_result$lengths)) {
      length_i <- rle_result$lengths[i]
      value_i <- rle_result$values[i]
      
      if (value_i && length_i >= 3) {  # 至少连续3个点
        significant_cleaned[current_pos:(current_pos + length_i - 1)] <- TRUE
      }
      current_pos <- current_pos + length_i
    }
    
    cat("    连续性清理后保留", sum(significant_cleaned, na.rm = TRUE), "个显著发展点\n")
    
    # 模型质量指标
    aic <- AIC(gam_model)
    deviance_explained <- model_summary$dev.expl
    r_squared <- model_summary$r.sq
    
    results <- list(
      age = age_range,
      predictions = predictions,
      pred_lower = pred_lower,
      pred_upper = pred_upper,
      derivatives_6mo_raw = derivatives_6mo,
      derivatives_annual = derivatives_annual,
      derivatives = derivatives_6mo_increment,  # 主要结果
      deriv_lower = deriv_lower,
      deriv_upper = deriv_upper,
      significant = significant_cleaned,  # 使用清理后的显著性
      significance_threshold = significance_threshold,  # 记录使用的阈值
      aic = aic,
      deviance_explained = deviance_explained,
      r_squared = r_squared,
      k_basis = k_basis,  # 修改：使用k_basis替代optimal_k
      residual_std = residual_std,
      n_bootstrap_actual = sum(valid_bootstrap),
      smooth_f_value = smooth_f_value,      
      smooth_p_value = smooth_p_value,      
      smooth_ref_df = smooth_ref_df,        
      smooth_edf = smooth_edf,              
      model = gam_model,
      model_summary = model_summary
    )
    
    return(results)
    
  }, error = function(e) {
    cat("    GAM fitting failed for", outcome_var, ":", e$message, "\n")
    return(NULL)
  })
}

# 创建自定义颜色映射
create_custom_colormap <- function() {
  colorRampPalette(c("#08306B", "#4878B6", "#FFFFFF", "#D73027", "#67000D"))(256)
}

# 绘制梯度栅格图
plot_gradient_raster_figure <- function(age_data, derivative_data_list, significance_list, metric_names) {
  
  # 准备数据
  all_vals <- unlist(derivative_data_list)
  abs_max <- max(abs(range(all_vals, na.rm = TRUE)))
  
  # 创建数据框
  raster_data <- data.frame()
  
  for (m_idx in seq_along(metric_names)) {
    for (a_idx in seq_along(age_data)) {
      if (significance_list[[m_idx]][a_idx]) {
        raster_data <- rbind(raster_data, data.frame(
          age = age_data[a_idx],
          metric = factor(metric_names[m_idx], levels = rev(metric_names)),
          value = derivative_data_list[[m_idx]][a_idx],
          significant = TRUE
        ))
      } else {
        raster_data <- rbind(raster_data, data.frame(
          age = age_data[a_idx],
          metric = factor(metric_names[m_idx], levels = rev(metric_names)),
          value = 0,
          significant = FALSE
        ))
      }
    }
  }
  
  # 创建图形
  p <- ggplot(raster_data) +
    geom_tile(aes(x = age, y = metric, fill = ifelse(significant, value, NA)), 
              width = 0.1, height = 0.8) +
    geom_tile(data = raster_data[!raster_data$significant, ],
              aes(x = age, y = metric), 
              fill = NA, color = "grey", linetype = "dashed", size = 0.2,
              width = 0.1, height = 0.8) +
    scale_fill_gradient2(low = "#08306B", mid = "#FFFFFF", high = "#D73027",
                         midpoint = 0, na.value = "transparent",
                         name = "6-month developmental\nincrement (sd units)",
                         limits = c(-abs_max, abs_max)) +
    scale_x_continuous(breaks = 6:19, limits = c(6, 19)) +
    scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    labs(x = "Age (years)", y = "") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text.y = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(face = "bold", size = 11),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 9)
    )
  
  return(p)
}

# 移除异常值函数
remove_outliers_4sd <- function(data, column) {
  mean_val <- mean(data[[column]], na.rm = TRUE)
  sd_val <- sd(data[[column]], na.rm = TRUE)
  lower_bound <- mean_val - 4 * sd_val
  upper_bound <- mean_val + 4 * sd_val
  
  return(data[data[[column]] >= lower_bound & data[[column]] <= upper_bound, ])
}

# 结果打印函数（修改后版本）
print_comprehensive_results <- function(gam_results, inverted_u_stats,
                                        stage_stats, completion_stats,
                                        successful_metrics) {
  
  cat("\n=== 修改后的显著性检验结果 ===\n")
  cat("显著性标准: |derivatives_6mo_increment| >=", significance_threshold, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  for (metric in successful_metrics) {
    results <- gam_results[[metric]]
    cat("\n", metric, ":\n", sep = "")
    cat("  AIC:", round(results$aic, 2), "\n")
    cat("  R²:", round(results$r_squared * 100, 1), "%\n")
    cat("  EDF:", round(results$smooth_edf, 2), "\n")
    cat("  F值:", round(results$smooth_f_value, 1), "\n")
    cat("  p值:", round(results$smooth_p_value, 4), "\n")
    cat("  显著发展点数量:", sum(results$significant, na.rm = TRUE), "\n")
    cat("  显著发展占比:", round(sum(results$significant, na.rm = TRUE) / length(results$significant) * 100, 1), "%\n")
    
    # 判断显著性水平
    if (results$smooth_p_value < 0.001) {
      cat("  平滑项显著性: *** (p < 0.001)\n")
    } else if (results$smooth_p_value < 0.01) {
      cat("  平滑项显著性: ** (p < 0.01)\n")
    } else if (results$smooth_p_value < 0.05) {
      cat("  平滑项显著性: * (p < 0.05)\n")
    } else if (results$smooth_p_value < 0.1) {
      cat("  平滑项显著性: . (p < 0.1)\n")
    } else {
      cat("  平滑项显著性: 不显著 (p >= 0.1)\n")
    }
  }
  
  cat("\n=== 倒U型发展模式分析 ===\n")
  cat(rep("=", 60), "\n", sep = "")
  
  inverted_u_metrics <- c()
  
  for (metric in successful_metrics) {
    stats <- inverted_u_stats[[metric]]
    if (stats$has_inverted_u) {
      cat("\n", metric, " - 检测到倒U型发展模式:\n", sep = "")
      cat("  峰值年龄: ", round(stats$peak_age, 1), "岁\n", sep = "")
      cat("  峰值水平: ", round(stats$peak_value, 3), " (Z分数)\n", sep = "")
      
      if (!any(is.na(stats$increase_period))) {
        cat("  上升期: ", round(stats$increase_period[1], 1), " - ",
            round(stats$increase_period[2], 1), "岁 (持续",
            round(stats$increase_duration, 1), "年)\n", sep = "")
      }
      
      if (!any(is.na(stats$decrease_period))) {
        cat("  下降期: ", round(stats$decrease_period[1], 1), " - ",
            round(stats$decrease_period[2], 1), "岁 (持续",
            round(stats$decrease_duration, 1), "年)\n", sep = "")
      }
      
      inverted_u_metrics <- c(inverted_u_metrics, metric)
    } else {
      cat("\n", metric, " - 未检测到明显的倒U型模式\n", sep = "")
    }
  }
  
  # 发展阶段变化率分析（修改：增加d_value作为核心指标）
  cat("\n=== 发展阶段变化率比较分析（核心执行功能指标）===\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # 修改核心抑制指标，增加d_value
  core_inhibition_metrics <- successful_metrics[
    grepl("Stroop_Incongruent|GoNoGo_Nogo_ACC|GoNoGo_d_value", successful_metrics)
  ]
  
  if (length(core_inhibition_metrics) > 0) {
    cat("分析的核心执行功能指标: ", paste(core_inhibition_metrics, collapse = ", "), "\n")
    
    for (metric in core_inhibition_metrics) {
      stats <- stage_stats[[metric]]
      cat("\n", metric, ":\n", sep = "")
      
      cat("  早期快速发展阶段 (6-8岁):\n")
      cat("    年龄范围: ", stats$early_rapid_6_8_age_range, "\n")
      cat("    平均变化率: ", round(stats$early_rapid_6_8_mean_rate, 3), " 标准差单位/年\n")
      cat("    平均绝对变化率: ", round(stats$early_rapid_6_8_abs_mean_rate, 3), " 标准差单位/年\n")
      cat("    持续时间: ", round(stats$early_rapid_6_8_duration, 1), "年 (",
          round(stats$early_rapid_6_8_duration_pct, 1), "%)\n")
      
      cat("  中期持续发展阶段 (9-12岁):\n")
      cat("    年龄范围: ", stats$middle_sustained_9_12_age_range, "\n")
      cat("    平均变化率: ", round(stats$middle_sustained_9_12_mean_rate, 3), " 标准差单位/年\n")
      cat("    平均绝对变化率: ", round(stats$middle_sustained_9_12_abs_mean_rate, 3), " 标准差单位/年\n")
      cat("    持续时间: ", round(stats$middle_sustained_9_12_duration, 1), "年 (",
          round(stats$middle_sustained_9_12_duration_pct, 1), "%)\n")
      
      cat("  后期发展阶段 (13-15岁):\n")
      cat("    年龄范围: ", stats$late_plateau_13_15_age_range, "\n")
      cat("    平均变化率: ", round(stats$late_plateau_13_15_mean_rate, 3), " 标准差单位/年\n")
      cat("    平均绝对变化率: ", round(stats$late_plateau_13_15_abs_mean_rate, 3), " 标准差单位/年\n")
      cat("    持续时间: ", round(stats$late_plateau_13_15_duration, 1), "年 (",
          round(stats$late_plateau_13_15_duration_pct, 1), "%)\n")
      
      cat("  稳定阶段 (16-18岁):\n")
      cat("    年龄范围: ", stats$final_stable_16_18_age_range, "\n")
      cat("    平均变化率: ", round(stats$final_stable_16_18_mean_rate, 3), " 标准差单位/年\n")
      cat("    平均绝对变化率: ", round(stats$final_stable_16_18_abs_mean_rate, 3), " 标准差单位/年\n")
      cat("    持续时间: ", round(stats$final_stable_16_18_duration, 1), "年 (",
          round(stats$final_stable_16_18_duration_pct, 1), "%)\n")
    }
  } else {
    cat("未找到核心执行功能指标\n")
  }
  
  # 多年龄切点发展完成度分析
  cat("\n=== 多年龄切点发展完成度分析 ===\n")
  cat(rep("=", 60), "\n", sep = "")
  
  cutoff_ages <- c(12, 13, 14, 15, 16)
  
  for (metric in successful_metrics) {
    stats <- completion_stats[[metric]]
    cat("\n", metric, ":\n", sep = "")
    
    if (!is.na(stats$sig_start)) {
      cat("  显著发展期: ", round(stats$sig_start, 1), " - ",
          round(stats$sig_end, 1), "岁\n", sep = "")
      cat("  显著发展总时长: ", round(stats$sig_duration, 1), "年\n", sep = "")
      
      for (cutoff in cutoff_ages) {
        completion_key <- paste0("completion_", cutoff)
        duration_key <- paste0("pre_", cutoff, "_duration")
        if (completion_key %in% names(stats)) {
          cat("  ", cutoff, "岁前完成度: ", round(stats[[completion_key]], 1),
              "% (发展时长: ", round(stats[[duration_key]], 1), "年)\n", sep = "")
        }
      }
    } else {
      cat("  未检测到显著发展期\n")
    }
  }
}

# ================================
# 主分析
# ================================

cat("=== 加载纵向数据文件 ===\n")

library(readxl)

# 在 Windows 下，路径可以用双反斜杠 \\ 或者正斜杠 /
tryCatch({
  df_stroop_baseline <- read_excel("C:/Users/Administrator/Desktop/3.Gonogo和Stroop比-主站点后测/Stroop-followup-matched.xlsx")
  df_stroop_followup <- read_excel("C:/Users/Administrator/Desktop/3.Gonogo和Stroop比-主站点后测/Stroop-followup-6mo.xlsx")
  df_gonogo_baseline <- read_excel("C:/Users/Administrator/Desktop/3.Gonogo和Stroop比-主站点后测/Gonogo-followup-matched.xlsx")
  df_gonogo_followup <- read_excel("C:/Users/Administrator/Desktop/3.Gonogo和Stroop比-主站点后测/Gonogo-followup-6mo.xlsx")
  cat("纵向数据文件加载成功\n")
}, error = function(e) {
  cat("加载文件时出错：\n", e$message, "\n")
})

# ================================
# 数据采样处理
# ================================

if (USE_SAMPLE_DATA) {
  cat("\n=== 使用", SAMPLE_PERCENTAGE * 100, "%的数据进行快速测试 ===\n")
  
  # 设置随机种子以确保可重现性
  set.seed(123)
  
  # 对每个数据集进行采样
  n_stroop_baseline <- max(1, floor(nrow(df_stroop_baseline) * SAMPLE_PERCENTAGE))
  n_stroop_followup <- max(1, floor(nrow(df_stroop_followup) * SAMPLE_PERCENTAGE))
  n_gonogo_baseline <- max(1, floor(nrow(df_gonogo_baseline) * SAMPLE_PERCENTAGE))
  n_gonogo_followup <- max(1, floor(nrow(df_gonogo_followup) * SAMPLE_PERCENTAGE))
  
  df_stroop_baseline <- df_stroop_baseline[sample(nrow(df_stroop_baseline), n_stroop_baseline), ]
  df_stroop_followup <- df_stroop_followup[sample(nrow(df_stroop_followup), n_stroop_followup), ]
  df_gonogo_baseline <- df_gonogo_baseline[sample(nrow(df_gonogo_baseline), n_gonogo_baseline), ]
  df_gonogo_followup <- df_gonogo_followup[sample(nrow(df_gonogo_followup), n_gonogo_followup), ]
  
  cat("采样后数据量:\n")
  cat("Stroop baseline:", nrow(df_stroop_baseline), "条观察\n")
  cat("Stroop followup:", nrow(df_stroop_followup), "条观察\n")
  cat("GoNoGo baseline:", nrow(df_gonogo_baseline), "条观察\n")
  cat("GoNoGo followup:", nrow(df_gonogo_followup), "条观察\n")
} else {
  cat("\n=== 使用完整数据集 ===\n")
}

# 添加时间点标签
df_stroop_baseline$timepoint <- "baseline"
df_stroop_followup$timepoint <- "followup"
df_gonogo_baseline$timepoint <- "baseline"
df_gonogo_followup$timepoint <- "followup"

# 合并数据
df_stroop_long <- rbind(df_stroop_baseline, df_stroop_followup)
df_gonogo_long <- rbind(df_gonogo_baseline, df_gonogo_followup)

cat("Stroop纵向数据:", nrow(df_stroop_long), "条观察\n")
cat("GoNoGo纵向数据:", nrow(df_gonogo_long), "条观察\n")

# 定义指标（修改：增加d_value）
stroop_metrics <- c("Incongruent", "Congruent", "Neutral", "Interference")
available_stroop <- intersect(stroop_metrics, names(df_stroop_long))

# 修改：增加d_value作为GoNoGo的核心指标
gonogo_metrics <- c("Go_ACC", "Nogo_ACC", "Go_RT_ms", "d_value")
available_gonogo <- intersect(gonogo_metrics, names(df_gonogo_long))

cat("使用Stroop指标:", paste(available_stroop, collapse = ", "), "\n")
cat("使用GoNoGo指标:", paste(available_gonogo, collapse = ", "), "\n")

# 创建综合纵向数据集
cat("\n=== 创建综合纵向数据集 ===\n")

all_data <- list()

# 处理Stroop数据
for (metric in available_stroop) {
  if (metric %in% names(df_stroop_long)) {
    temp_df <- df_stroop_long[, c("Subject", "Age", metric, "timepoint")]
    names(temp_df)[names(temp_df) == metric] <- "value"
    temp_df$task <- "Stroop"
    temp_df$metric <- metric
    temp_df$task_metric <- paste0("Stroop_", metric)
    all_data[[length(all_data) + 1]] <- temp_df
  }
}

# 处理GoNoGo数据（修改：包含d_value）
for (metric in available_gonogo) {
  if (metric %in% names(df_gonogo_long)) {
    temp_df <- df_gonogo_long[, c("Subject", "Age", metric, "timepoint")]
    names(temp_df)[names(temp_df) == metric] <- "value"
    temp_df$task <- "GoNoGo"
    temp_df$metric <- metric
    temp_df$task_metric <- paste0("GoNoGo_", metric)
    all_data[[length(all_data) + 1]] <- temp_df
  }
}

# 合并所有数据
df_combined <- do.call(rbind, all_data)
df_combined$age <- as.numeric(df_combined$Age)

# 移除缺失值
df_combined <- df_combined[!is.na(df_combined$value) & !is.na(df_combined$age), ]
cat("合并纵向数据集:", nrow(df_combined), "条观察，涵盖", 
    length(unique(df_combined$task_metric)), "个指标\n")

# Z分数标准化
df_combined <- df_combined %>%
  group_by(task_metric) %>%
  mutate(value_z = scale(value)[, 1]) %>%
  ungroup()

# 数据质量控制
cat("原始样本量:", nrow(df_combined), "\n")

# 移除异常值
df_clean <- df_combined
for (task_metric in unique(df_clean$task_metric)) {
  task_data <- df_clean[df_clean$task_metric == task_metric, ]
  task_data_clean <- remove_outliers_4sd(task_data, "value_z")
  df_clean <- df_clean[!(df_clean$task_metric == task_metric & 
                           !df_clean$Subject %in% task_data_clean$Subject), ]
}

cat("异常值移除后:", nrow(df_clean), "条观察\n")

# 定义年龄范围
age_range <- seq(min(df_clean$age, na.rm = TRUE), 
                 max(df_clean$age, na.rm = TRUE), 
                 by = 0.1)

# GAM分析（修改后方法）
cat("\n=== 所有指标的GAM分析（使用修改后显著性标准）===\n")


gam_results <- list()
successful_metrics <- character(0)


for (task_metric in unique(df_clean$task_metric)) {
  cat("\n分析", task_metric, "...\n")
  metric_data <- df_clean[df_clean$task_metric == task_metric, ]
  
  if (nrow(metric_data) < 30) {
    cat("  跳过", task_metric, ": 数据不足 (", nrow(metric_data), "条观察)\n")
    next
  }
  
  # 检查时间点
  timepoints <- unique(metric_data$timepoint)
  if (length(timepoints) < 2) {
    cat("  跳过", task_metric, ": 只有一个时间点可用\n")
    next
  }
  
  # 确保Subject是因子
  metric_data$Subject <- as.factor(metric_data$Subject)
  metric_data$timepoint <- as.factor(metric_data$timepoint)
  
  results <- compute_gam_longitudinal(
    metric_data, "value_z", "age", "Subject", "timepoint",
    age_range, n_bootstrap = 3000, k_basis = 10, 
    significance_threshold = significance_threshold  # 传递新的阈值
  )
  
  if (!is.null(results)) {
    gam_results[[task_metric]] <- results
    successful_metrics <- c(successful_metrics, task_metric)
    cat("  结果: AIC:", round(results$aic, 2), ", R²:", round(results$r_squared * 100, 1), 
        "%, k:", results$k_basis, ", F值:", round(results$smooth_f_value, 1),
        ", p值:", round(results$smooth_p_value, 4), "\n")
    cat("  显著发展点:", sum(results$significant, na.rm = TRUE), "/", length(results$significant), 
        " (", round(sum(results$significant, na.rm = TRUE) / length(results$significant) * 100, 1), "%)\n")
  }
}

if (length(successful_metrics) == 0) {
  cat("没有成功的GAM分析。退出。\n")
  stop()
}

cat("\n成功分析:", length(successful_metrics), "个指标\n")

# ================================
# 发展模式分析
# ================================

cat("\n=== 检测倒U型发展模式 ===\n")
inverted_u_stats <- list()

for (task_metric in successful_metrics) {
  results <- gam_results[[task_metric]]
  inverted_u_stats[[task_metric]] <- detect_inverted_u_pattern(results, age_range)
}

cat("\n=== 计算发展阶段变化率 ===\n")
stage_stats <- list()

for (task_metric in successful_metrics) {
  results <- gam_results[[task_metric]]
  stage_stats[[task_metric]] <- calculate_development_stage_rates(results, age_range)
}

cat("\n=== 计算发展完成度 ===\n")
completion_stats <- list()
cutoff_ages <- c(12, 13, 14, 15, 16)

for (task_metric in successful_metrics) {
  results <- gam_results[[task_metric]]
  completion_stats[[task_metric]] <- calculate_development_completion(results, age_range, cutoff_ages)
}

# 打印综合结果（使用修改后的函数）
print_comprehensive_results(gam_results, inverted_u_stats, stage_stats,
                            completion_stats, successful_metrics, significance_threshold)

# ================================
# 可视化
# ================================

cat("\n=== 创建可视化 ===\n")

# 1. 轨迹图
cat("创建发育轨迹图...\n")

# 准备轨迹数据
trajectory_data <- data.frame()
for (metric in successful_metrics) {
  results <- gam_results[[metric]]
  temp_data <- data.frame(
    age = results$age,
    prediction = results$predictions,
    pred_lower = results$pred_lower,
    pred_upper = results$pred_upper,
    task_metric = metric,
    task = ifelse(grepl("Stroop", metric), "Stroop", "GoNoGo"),
    k_basis = results$k_basis  # 修改：使用k_basis
  )
  trajectory_data <- rbind(trajectory_data, temp_data)
}

# 绘制轨迹图
trajectory_plot <- ggplot(trajectory_data, aes(x = age, y = prediction)) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = task), alpha = 0.2) +
  geom_line(aes(color = task, linetype = task_metric), linewidth = 1.2) +
  scale_color_manual(values = c("Stroop" = "#1f77b4", "GoNoGo" = "#ff7f0e")) +
  scale_fill_manual(values = c("Stroop" = "#1f77b4", "GoNoGo" = "#ff7f0e")) +
  scale_x_continuous(breaks = 6:19, limits = c(6, 19)) +
  labs(x = "Age (years)", y = "Z-scored Performance",
       color = "Task", fill = "Task", linetype = "Metric") +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUTPUT_PATH, "developmental_trajectories_modified.png"), trajectory_plot, width = 6, height = 4, dpi = 800)

# 2. 梯度栅格图
if (length(successful_metrics) > 1) {
  cat("创建6个月发育增量栅格图（修改后显著性标准）...\n")
  
  derivative_data_list <- lapply(successful_metrics, function(m) gam_results[[m]]$derivatives)
  significance_list <- lapply(successful_metrics, function(m) gam_results[[m]]$significant)
  
  raster_plot <- plot_gradient_raster_figure(age_range, derivative_data_list, 
                                             significance_list, successful_metrics)
  
  # 添加标注显示使用的新标准
  raster_plot <- raster_plot + 
    labs(subtitle = paste("Significance threshold: |increment| ≥", significance_threshold))
  
  ggsave(file.path(OUTPUT_PATH, "developmental_increment_raster_modified.png"), raster_plot, width = 7, height = 3.5, dpi = 800)
}

# ================================
# Enhanced Summary Statistics Section (修改后版本)
# ================================

cat("\n=== 汇总统计（修改后版本） ===\n")

# 创建增强汇总表
summary_stats <- data.frame()

for (task_metric in successful_metrics) {
  results <- gam_results[[task_metric]]
  sig_periods <- results$significant
  
  # 基本统计
  stats_row <- data.frame(
    Metric = task_metric,
    Task = strsplit(task_metric, "_")[[1]][1],
    Variable = paste(strsplit(task_metric, "_")[[1]][-1], collapse = "_"),
    N_Observations = nrow(df_clean[df_clean$task_metric == task_metric, ]),
    AIC = results$aic,
    R_squared = results$r_squared,
    Deviance_Explained = results$deviance_explained,
    K_Basis = results$k_basis,  # 修改：使用k_basis
    Residual_Std = results$residual_std,
    Bootstrap_Samples = results$n_bootstrap_actual,
    # 平滑项显著性
    Smooth_F_Value = results$smooth_f_value,
    Smooth_P_Value = results$smooth_p_value,
    Smooth_Ref_DF = results$smooth_ref_df,
    Smooth_EDF = results$smooth_edf,
    # 修改后显著性统计
    Significance_Threshold = results$significance_threshold,
    Significant_Periods_Count = sum(sig_periods, na.rm = TRUE),
    Significant_Periods_Pct = sum(sig_periods, na.rm = TRUE) / length(sig_periods) * 100,
    Max_6mo_Increment = max(abs(results$derivatives), na.rm = TRUE),
    Mean_6mo_Increment = mean(results$derivatives, na.rm = TRUE),
    Max_Annual_Rate = max(abs(results$derivatives_annual), na.rm = TRUE),
    Mean_Annual_Rate = mean(results$derivatives_annual, na.rm = TRUE)
  )
  
  # 添加显著期信息
  if (any(sig_periods, na.rm = TRUE)) {
    sig_ages <- results$age[sig_periods]
    stats_row$Sig_Age_Start <- min(sig_ages, na.rm = TRUE)
    stats_row$Sig_Age_End <- max(sig_ages, na.rm = TRUE)
    stats_row$Sig_Duration <- max(sig_ages, na.rm = TRUE) - min(sig_ages, na.rm = TRUE)
  } else {
    stats_row$Sig_Age_Start <- NA
    stats_row$Sig_Age_End <- NA
    stats_row$Sig_Duration <- 0
  }
  
  # 峰值变化年龄
  max_change_idx <- which.max(abs(results$derivatives))
  stats_row$Peak_Change_Age <- results$age[max_change_idx]
  stats_row$Peak_Change_Value <- results$derivatives[max_change_idx]
  
  # 倒U型统计
  u_stats <- inverted_u_stats[[task_metric]]
  stats_row$Has_Inverted_U <- u_stats$has_inverted_u
  stats_row$U_Peak_Age <- ifelse(is.na(u_stats$peak_age), NA, u_stats$peak_age)
  stats_row$U_Peak_Value <- ifelse(is.na(u_stats$peak_value), NA, u_stats$peak_value)
  stats_row$U_Increase_Duration <- u_stats$increase_duration
  stats_row$U_Decrease_Duration <- u_stats$decrease_duration
  
  # 发展阶段统计
  stage_info <- stage_stats[[task_metric]]
  stats_row$Early_Stage_6_8_Rate <- stage_info$early_rapid_6_8_mean_rate
  stats_row$Early_Stage_6_8_Abs_Rate <- stage_info$early_rapid_6_8_abs_mean_rate
  stats_row$Middle_Stage_9_12_Rate <- stage_info$middle_sustained_9_12_mean_rate
  stats_row$Middle_Stage_9_12_Abs_Rate <- stage_info$middle_sustained_9_12_abs_mean_rate
  stats_row$Late_Stage_13_15_Rate <- stage_info$late_plateau_13_15_mean_rate
  stats_row$Late_Stage_13_15_Abs_Rate <- stage_info$late_plateau_13_15_abs_mean_rate
  stats_row$Final_Stage_16_18_Rate <- stage_info$final_stable_16_18_mean_rate
  stats_row$Final_Stage_16_18_Abs_Rate <- stage_info$final_stable_16_18_abs_mean_rate
  
  # 完成度统计
  completion_info <- completion_stats[[task_metric]]
  cutoff_ages <- c(12, 13, 14, 15, 16)
  for (cutoff in cutoff_ages) {
    completion_key <- paste0("completion_", cutoff)
    duration_key <- paste0("pre_", cutoff, "_duration")
    if (completion_key %in% names(completion_info)) {
      stats_row[[paste0("Completion_", cutoff, "_Pct")]] <- completion_info[[completion_key]]
      stats_row[[paste0("Pre_", cutoff, "_Duration")]] <- completion_info[[duration_key]]
    }
  }
  
  summary_stats <- rbind(summary_stats, stats_row)
}

# 修复打印问题：只对数值列应用round函数
numeric_cols <- sapply(summary_stats, is.numeric)
summary_stats_display <- summary_stats
summary_stats_display[numeric_cols] <- round(summary_stats_display[numeric_cols], 3)

print(summary_stats_display)

# ================================
# 保存结果（修改后版本）
# ================================

cat("\n=== 保存结果 ===\n")

# 1. 保存汇总统计表
write.csv(summary_stats_display, file.path(OUTPUT_PATH, "modified_gam_summary_statistics.csv"), row.names = FALSE)
cat("保存修改后汇总统计表到", file.path(OUTPUT_PATH, "modified_gam_summary_statistics.csv"), "\n")

# 2. 保存详细的GAM结果
detailed_results <- list()
for (metric in successful_metrics) {
  results <- gam_results[[metric]]
  detailed_results[[metric]] <- data.frame(
    age = results$age,
    predictions = results$predictions,
    pred_lower = results$pred_lower,
    pred_upper = results$pred_upper,
    derivatives_6mo_increment = results$derivatives,
    derivatives_annual = results$derivatives_annual,
    deriv_lower = results$deriv_lower,
    deriv_upper = results$deriv_upper,
    significant_original_method = (results$deriv_lower > 0) | (results$deriv_upper < 0),  # 原始方法
    significant_modified_method = results$significant,  # 修改后方法
    significance_threshold = results$significance_threshold
  )
}

# 3. 保存发展阶段分析结果
stage_analysis_results <- data.frame()
for (task_metric in successful_metrics) {
  stage_info <- stage_stats[[task_metric]]
  row_data <- data.frame(
    Metric = task_metric,
    Task = strsplit(task_metric, "_")[[1]][1],
    Variable = paste(strsplit(task_metric, "_")[[1]][-1], collapse = "_")
  )
  for (name in names(stage_info)) {
    row_data[[name]] <- stage_info[[name]]
  }
  stage_analysis_results <- rbind(stage_analysis_results, row_data)
}

# 4. 保存发展完成度分析结果
completion_analysis_results <- data.frame()
for (task_metric in successful_metrics) {
  completion_info <- completion_stats[[task_metric]]
  row_data <- data.frame(
    Metric = task_metric,
    Task = strsplit(task_metric, "_")[[1]][1],
    Variable = paste(strsplit(task_metric, "_")[[1]][-1], collapse = "_")
  )
  for (name in names(completion_info)) {
    row_data[[name]] <- completion_info[[name]]
  }
  completion_analysis_results <- rbind(completion_analysis_results, row_data)
}

# 5. 保存倒U型分析结果
inverted_u_results <- data.frame()
for (task_metric in successful_metrics) {
  u_info <- inverted_u_stats[[task_metric]]
  row_data <- data.frame(
    Metric = task_metric,
    Task = strsplit(task_metric, "_")[[1]][1],
    Variable = paste(strsplit(task_metric, "_")[[1]][-1], collapse = "_"),
    Has_Inverted_U = u_info$has_inverted_u,
    Peak_Age = ifelse(is.na(u_info$peak_age), NA, u_info$peak_age),
    Peak_Value = ifelse(is.na(u_info$peak_value), NA, u_info$peak_value),
    Increase_Duration = u_info$increase_duration,
    Decrease_Duration = u_info$decrease_duration
  )
  inverted_u_results <- rbind(inverted_u_results, row_data)
}

# 保存为Excel文件
library(openxlsx)
wb <- createWorkbook()

# 添加工作表
addWorksheet(wb, "Modified_Summary")
addWorksheet(wb, "Stage_Analysis")
addWorksheet(wb, "Completion_Analysis")
addWorksheet(wb, "Inverted_U_Analysis")
addWorksheet(wb, "Clean_Data")
addWorksheet(wb, "Method_Comparison")

# 写入数据
writeData(wb, "Modified_Summary", summary_stats_display)
writeData(wb, "Stage_Analysis", stage_analysis_results)
writeData(wb, "Completion_Analysis", completion_analysis_results)
writeData(wb, "Inverted_U_Analysis", inverted_u_results)
writeData(wb, "Clean_Data", df_clean)

# 创建方法比较表
method_comparison <- data.frame()
for (metric in successful_metrics) {
  results <- gam_results[[metric]]
  original_significant <- (results$deriv_lower > 0) | (results$deriv_upper < 0)
  modified_significant <- results$significant
  
  comparison_row <- data.frame(
    Metric = metric,
    Original_Method_Count = sum(original_significant, na.rm = TRUE),
    Original_Method_Pct = round(sum(original_significant, na.rm = TRUE) / length(original_significant) * 100, 2),
    Modified_Method_Count = sum(modified_significant, na.rm = TRUE),
    Modified_Method_Pct = round(sum(modified_significant, na.rm = TRUE) / length(modified_significant) * 100, 2),
    Difference_Count = sum(modified_significant, na.rm = TRUE) - sum(original_significant, na.rm = TRUE),
    Threshold_Used = results$significance_threshold
  )
  method_comparison <- rbind(method_comparison, comparison_row)
}

writeData(wb, "Method_Comparison", method_comparison)

# 为每个指标添加详细结果
for (metric in names(detailed_results)) {
  sheet_name <- substr(gsub("[^A-Za-z0-9]", "_", metric), 1, 31)  # Excel工作表名称限制
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, detailed_results[[metric]])
}

saveWorkbook(wb, file.path(OUTPUT_PATH, "modified_gam_detailed_results.xlsx"), overwrite = TRUE)
cat("保存修改后详细结果到", file.path(OUTPUT_PATH, "modified_gam_detailed_results.xlsx"), "\n")

# 6. 保存平滑项显著性汇总表（修改后版本）
significance_summary <- data.frame()
for (metric in successful_metrics) {
  results <- gam_results[[metric]]
  
  # 确定显著性等级
  p_val <- results$smooth_p_value
  if (is.na(p_val)) {
    sig_level <- ""
    sig_text <- "无数据"
  } else if (p_val < 0.001) {
    sig_level <- "***"
    sig_text <- "p < 0.001"
  } else if (p_val < 0.01) {
    sig_level <- "**"
    sig_text <- "p < 0.01"
  } else if (p_val < 0.05) {
    sig_level <- "*"
    sig_text <- "p < 0.05"
  } else if (p_val < 0.1) {
    sig_level <- "."
    sig_text <- "p < 0.1"
  } else {
    sig_level <- ""
    sig_text <- "不显著"
  }
  
  row_data <- data.frame(
    Metric = metric,
    Task = strsplit(metric, "_")[[1]][1],
    Variable = paste(strsplit(metric, "_")[[1]][-1], collapse = "_"),
    F_Value = round(results$smooth_f_value, 3),
    P_Value = round(results$smooth_p_value, 6),
    Ref_DF = round(results$smooth_ref_df, 2),
    EDF = round(results$smooth_edf, 2),
    Significance_Level = sig_level,
    Significance_Text = sig_text,
    R_squared = round(results$r_squared * 100, 1),
    AIC = round(results$aic, 2),
    Modified_Sig_Count = sum(results$significant, na.rm = TRUE),
    Modified_Sig_Pct = round(sum(results$significant, na.rm = TRUE) / length(results$significant) * 100, 1),
    Threshold_Used = results$significance_threshold,
    stringsAsFactors = FALSE
  )
  
  significance_summary <- rbind(significance_summary, row_data)
}

# 按任务和显著性排序
significance_summary <- significance_summary[order(significance_summary$Task, significance_summary$P_Value), ]

write.csv(significance_summary, file.path(OUTPUT_PATH, "modified_smooth_significance_summary.csv"), 
          row.names = FALSE, fileEncoding = "UTF-8")
cat("保存修改后平滑项显著性汇总表到", file.path(OUTPUT_PATH, "modified_smooth_significance_summary.csv"), "\n")

# 7. 保存模型质量诊断（修改后版本）
model_diagnostics <- data.frame()
for (metric in successful_metrics) {
  results <- gam_results[[metric]]
  
  # 基本模型信息
  diag_row <- data.frame(
    Metric = metric,
    Task = strsplit(metric, "_")[[1]][1],
    Variable = paste(strsplit(metric, "_")[[1]][-1], collapse = "_"),
    AIC = results$aic,
    R_squared = results$r_squared,
    Deviance_Explained = results$deviance_explained,
    Smooth_F_Value = results$smooth_f_value,
    Smooth_P_Value = results$smooth_p_value,
    Smooth_EDF = results$smooth_edf,
    K_Basis = results$k_basis,  # 修改：使用k_basis
    Residual_Std = results$residual_std,
    Bootstrap_Success_Rate = results$n_bootstrap_actual / 3000 * 100,  # 基于30次bootstrap
    Mean_Abs_Derivative = mean(abs(results$derivatives), na.rm = TRUE),
    Max_Abs_Derivative = max(abs(results$derivatives), na.rm = TRUE),
    Significance_Threshold = results$significance_threshold,
    Modified_Significant_Count = sum(results$significant, na.rm = TRUE),
    Modified_Significant_Age_Span = ifelse(any(results$significant, na.rm = TRUE), 
                                           max(results$age[results$significant], na.rm = TRUE) - 
                                             min(results$age[results$significant], na.rm = TRUE), 0),
    Modified_Significant_Periods_Pct = sum(results$significant, na.rm = TRUE) / length(results$significant) * 100
  )
  
  model_diagnostics <- rbind(model_diagnostics, diag_row)
}

write.csv(model_diagnostics, file.path(OUTPUT_PATH, "modified_model_diagnostics.csv"), row.names = FALSE)
cat("保存修改后模型诊断信息到", file.path(OUTPUT_PATH, "modified_model_diagnostics.csv"), "\n")

# 8. 创建年龄特异性变化率表（修改后版本）
age_specific_changes <- data.frame(age = age_range)
for (metric in successful_metrics) {
  results <- gam_results[[metric]]
  age_specific_changes[[paste0(metric, "_increment")]] <- results$derivatives
  age_specific_changes[[paste0(metric, "_original_significant")]] <- (results$deriv_lower > 0) | (results$deriv_upper < 0)
  age_specific_changes[[paste0(metric, "_modified_significant")]] <- results$significant
  age_specific_changes[[paste0(metric, "_abs_increment")]] <- abs(results$derivatives)
  age_specific_changes[[paste0(metric, "_meets_threshold")]] <- abs(results$derivatives) >= results$significance_threshold
}

write.csv(age_specific_changes, file.path(OUTPUT_PATH, "modified_age_specific_developmental_changes.csv"), row.names = FALSE)
cat("保存修改后年龄特异性变化表到", file.path(OUTPUT_PATH, "modified_age_specific_developmental_changes.csv"), "\n")

# ================================
# 增强的数据可视化（修改后版本）
# ================================

cat("\n=== 创建增强可视化（修改后版本） ===\n")

# 修复ggplot2警告：使用linewidth替代size
# 更新轨迹图
trajectory_plot_fixed <- ggplot(trajectory_data, aes(x = age, y = prediction)) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = task), alpha = 0.2) +
  geom_line(aes(color = task, linetype = task_metric), linewidth = 1.2) +
  scale_color_manual(values = c("Stroop" = "#1f77b4", "GoNoGo" = "#ff7f0e")) +
  scale_fill_manual(values = c("Stroop" = "#1f77b4", "GoNoGo" = "#ff7f0e")) +
  scale_x_continuous(breaks = 6:19, limits = c(6, 19)) +
  labs(x = "Age (years)", y = "Z-scored Performance",
       color = "Task", fill = "Task", linetype = "Metric",
       title = "Developmental Trajectories of Cognitive Performance (Modified Significance)",
       subtitle = paste("Significance threshold: |increment| ≥", significance_threshold)) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

ggsave(file.path(OUTPUT_PATH, "modified_enhanced_developmental_trajectories.png"), trajectory_plot_fixed, 
       width = 7, height = 5, dpi = 800)

# 创建每个指标的单独详细图（修改后版本）
individual_plots_dir <- file.path(OUTPUT_PATH, "modified_individual_metric_plots")
if (!dir.exists(individual_plots_dir)) {
  dir.create(individual_plots_dir, recursive = TRUE)
}

for (metric in successful_metrics) {
  results <- gam_results[[metric]]
  
  # 预测轨迹图
  pred_data <- data.frame(
    age = results$age,
    prediction = results$predictions,
    pred_lower = results$pred_lower,
    pred_upper = results$pred_upper
  )
  
  pred_plot <- ggplot(pred_data, aes(x = age)) +
    geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.15, fill = "blue") +
    geom_line(aes(y = prediction), color = "blue", linewidth = 1.5) +
    scale_x_continuous(breaks = 6:19, limits = c(6, 19)) +
    labs(x = "Age (years)", y = "Z-scored Performance",
         title = paste("Trajectory:", gsub("_", " ", metric)),
         subtitle = paste0("AIC: ", round(results$aic, 1), ", R²: ", 
                           round(results$r_squared * 100, 1), "%, F: ", 
                           round(results$smooth_f_value, 1), ", p: ", 
                           round(results$smooth_p_value, 4))) +
    theme_minimal() +
    theme(
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  
  # 发育增量图（修改后版本）
  deriv_data <- data.frame(
    age = results$age,
    increment = results$derivatives,
    increment_lower = results$deriv_lower,
    increment_upper = results$deriv_upper,
    original_significant = (results$deriv_lower > 0) | (results$deriv_upper < 0),
    modified_significant = results$significant,
    abs_increment = abs(results$derivatives),
    threshold = results$significance_threshold
  )
  
  deriv_plot <- ggplot(deriv_data, aes(x = age)) +
    geom_ribbon(aes(ymin = increment_lower, ymax = increment_upper), 
                alpha = 0.15, fill = "#006d2c") +
    geom_line(aes(y = increment), color = "#006d2c", linewidth = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = c(results$significance_threshold, -results$significance_threshold), 
               linetype = "dotted", color = "red", alpha = 0.7) +
    geom_point(data = deriv_data[deriv_data$modified_significant, ], 
               aes(y = increment), color = "red", size = 1.5) +
    geom_point(data = deriv_data[deriv_data$original_significant & !deriv_data$modified_significant, ], 
               aes(y = increment), color = "orange", size = 1.2, alpha = 0.6) +
    scale_x_continuous(breaks = 6:19, limits = c(6, 19)) +
    labs(x = "Age (years)", y = "6-Month Developmental\nIncrement (SD units)",
         title = paste("Growth Rate:", gsub("_", " ", metric)),
         subtitle = paste("Red dots: modified significant (|increment| ≥", results$significance_threshold, 
                          "); Orange: original method only; Red lines: threshold")) +
    theme_minimal() +
    theme(
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  
  # 组合图
  combined_plot <- gridExtra::grid.arrange(pred_plot, deriv_plot, ncol = 1)
  
  filename <- file.path(individual_plots_dir, paste0(gsub("[^A-Za-z0-9]", "_", metric), "_modified_detailed.png"))
  ggsave(filename, combined_plot, width = 8, height = 8, dpi = 800)
}

cat("保存", length(successful_metrics), "个修改后详细指标图到", individual_plots_dir, "\n")

# 创建方法比较可视化
if (length(successful_metrics) >= 2) {
  
  # 准备方法比较数据
  comparison_data <- data.frame()
  
  for (metric in successful_metrics) {
    results <- gam_results[[metric]]
    original_significant <- (results$deriv_lower > 0) | (results$deriv_upper < 0)
    modified_significant <- results$significant
    
    for (i in 1:length(results$age)) {
      comparison_data <- rbind(comparison_data, data.frame(
        Metric = metric,
        Task = strsplit(metric, "_")[[1]][1],
        Age = results$age[i],
        Increment = results$derivatives[i],
        Abs_Increment = abs(results$derivatives[i]),
        Original_Significant = original_significant[i],
        Modified_Significant = modified_significant[i],
        Method_Agreement = original_significant[i] == modified_significant[i],
        Threshold = results$significance_threshold
      ))
    }
  }
  
  # 绘制方法比较散点图
  method_comparison_plot <- ggplot(comparison_data, aes(x = Age, y = Abs_Increment)) +
    geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed", linewidth = 1) +
    geom_point(aes(color = interaction(Original_Significant, Modified_Significant),
                   shape = Task), alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c(
        "FALSE.FALSE" = "lightgray",    # 两种方法都不显著
        "TRUE.FALSE" = "orange",        # 仅原始方法显著
        "FALSE.TRUE" = "blue",          # 仅修改方法显著
        "TRUE.TRUE" = "red"             # 两种方法都显著
      ),
      labels = c(
        "FALSE.FALSE" = "Both non-significant",
        "TRUE.FALSE" = "Original only",
        "FALSE.TRUE" = "Modified only", 
        "TRUE.TRUE" = "Both significant"
      ),
      name = "Significance Status"
    ) +
    scale_shape_manual(values = c("Stroop" = 16, "GoNoGo" = 17)) +
    scale_x_continuous(breaks = 6:19, limits = c(6, 19)) +
    facet_wrap(~ Metric, scales = "free_y") +
    labs(x = "Age (years)", y = "Absolute 6-Month Increment",
         title = "Comparison of Original vs Modified Significance Methods",
         subtitle = paste("Red dashed line: modified threshold (≥", significance_threshold, ")"),
         shape = "Task") +
    theme_minimal() +
    theme(
      axis.title = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      strip.text = element_text(size = 8),
      legend.position = "bottom"
    )
  
  ggsave(file.path(OUTPUT_PATH, "method_comparison_scatterplot.png"), method_comparison_plot, 
         width = 12, height = 8, dpi = 800)
}

# 创建汇总比较条形图
summary_comparison_data <- data.frame()
for (metric in successful_metrics) {
  results <- gam_results[[metric]]
  original_count <- sum((results$deriv_lower > 0) | (results$deriv_upper < 0), na.rm = TRUE)
  modified_count <- sum(results$significant, na.rm = TRUE)
  
  summary_comparison_data <- rbind(summary_comparison_data, data.frame(
    Metric = metric,
    Task = strsplit(metric, "_")[[1]][1],
    Original_Count = original_count,
    Modified_Count = modified_count,
    Original_Pct = round(original_count / length(results$significant) * 100, 1),
    Modified_Pct = round(modified_count / length(results$significant) * 100, 1),
    Difference = modified_count - original_count
  ))
}

# 重塑数据用于绘图
comparison_long <- reshape2::melt(summary_comparison_data[, c("Metric", "Task", "Original_Pct", "Modified_Pct")],
                                  id.vars = c("Metric", "Task"),
                                  variable.name = "Method",
                                  value.name = "Percentage")

comparison_bar_plot <- ggplot(comparison_long, aes(x = Metric, y = Percentage, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Original_Pct" = "#1f77b4", "Modified_Pct" = "#ff7f0e"),
                    labels = c("Original_Pct" = "Original Method", "Modified_Pct" = "Modified Method")) +
  labs(x = "Metric", y = "Percentage of Significant Time Points (%)",
       title = "Comparison of Significance Detection Methods",
       subtitle = paste("Modified threshold: |increment| ≥", significance_threshold),
       fill = "Method") +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_PATH, "method_comparison_barplot.png"), comparison_bar_plot, 
       width = 8, height = 6, dpi = 800)

# ================================
# 生成修改后分析报告
# ================================

cat("\n=== 生成修改后分析报告 ===\n")

# 创建markdown报告
report_content <- paste0(
  "# GAM纵向发育分析报告（修改后显著性标准，包含d_value）\n\n",
  "**分析日期**: ", Sys.Date(), "\n\n",
  ifelse(USE_SAMPLE_DATA, 
         paste0("**数据类型**: 使用", SAMPLE_PERCENTAGE * 100, "%采样数据进行快速测试\n\n"),
         "**数据类型**: 完整数据集\n\n"),
  "**分析方法**: 标准GAM方法，使用REML平滑参数选择，k=10基函数\n\n",
  "**重要修改**: 显著性判断标准改为 |derivatives_6mo_increment| ≥ ", significance_threshold, "\n\n",
  "**指标更新**: GoNoGo任务增加d'值(d_value)作为核心抑制指标\n\n",
  "**数据概要**: \n",
  "- 总观察数: ", nrow(df_clean), "\n",
  "- 成功分析指标数: ", length(successful_metrics), "\n",
  "- 年龄范围: ", round(min(age_range), 1), " - ", round(max(age_range), 1), " 岁\n",
  "- 显著性阈值: ", significance_threshold, " 标准差单位\n",
  "- 核心抑制指标: Stroop_Incongruent, GoNoGo_Nogo_ACC, GoNoGo_d_value\n\n",
  
  "## 主要发现\n\n",
  "### 修改后显著性检验结果\n"
)

# 添加每个指标的详细信息
for (i in 1:nrow(summary_stats_display)) {
  metric_info <- summary_stats_display[i, ]
  significance_level <- ifelse(metric_info$Smooth_P_Value < 0.001, "***",
                               ifelse(metric_info$Smooth_P_Value < 0.01, "**",
                                      ifelse(metric_info$Smooth_P_Value < 0.05, "*", "")))
  
  # 标记核心指标
  is_core_inhibition <- grepl("Stroop_Incongruent|GoNoGo_Nogo_ACC|GoNoGo_d_value", metric_info$Metric)
  core_marker <- ifelse(is_core_inhibition, " [核心抑制指标]", "")
  
  report_content <- paste0(
    report_content,
    "**", metric_info$Metric, "**", significance_level, core_marker, ":\n",
    "- AIC: ", metric_info$AIC, ", R²: ", metric_info$R_squared, "%, EDF: ", metric_info$Smooth_EDF, "\n",
    "- F值: ", metric_info$Smooth_F_Value, ", p值: ", metric_info$Smooth_P_Value, "\n",
    "- 修改后显著发育点数量: ", metric_info$Significant_Periods_Count, "\n",
    "- 修改后显著发育期占比: ", metric_info$Significant_Periods_Pct, "%\n",
    "- 最大6个月增量: ", metric_info$Max_6mo_Increment, " SD\n",
    "- 峰值变化年龄: ", metric_info$Peak_Change_Age, " 岁\n\n"
  )
}

# 添加方法比较总结
report_content <- paste0(
  report_content,
  "\n### 方法比较总结\n",
  "修改后的显著性标准 (|increment| ≥ ", significance_threshold, ") 与原始方法比较：\n\n"
)

for (i in 1:nrow(method_comparison)) {
  row <- method_comparison[i, ]
  report_content <- paste0(
    report_content,
    "**", row$Metric, "**:\n",
    "- 原始方法显著点: ", row$Original_Method_Count, " (", row$Original_Method_Pct, "%)\n",
    "- 修改方法显著点: ", row$Modified_Method_Count, " (", row$Modified_Method_Pct, "%)\n",
    "- 差异: ", ifelse(row$Difference_Count >= 0, "+", ""), row$Difference_Count, " 个点\n\n"
  )
}

# 添加核心抑制指标特别分析
core_inhibition_metrics_found <- successful_metrics[
  grepl("Stroop_Incongruent|GoNoGo_Nogo_ACC|GoNoGo_d_value", successful_metrics)
]

if (length(core_inhibition_metrics_found) > 0) {
  report_content <- paste0(
    report_content,
    "\n### 核心抑制指标特别分析\n",
    "本次分析包含以下核心抑制控制指标：\n",
    "- **Stroop_Incongruent**: 冲突抑制能力\n",
    "- **GoNoGo_Nogo_ACC**: 反应抑制准确性\n",
    "- **GoNoGo_d_value**: 信号检测理论的辨别力指标\n\n",
    
    "d'值(d_value)作为新增的核心指标，反映了被试在Go/NoGo任务中区分Go信号和NoGo信号的能力，",
    "是比简单准确率更敏感的抑制控制测量指标。\n\n"
  )
}

# 添加发展阶段分析总结
report_content <- paste0(
  report_content,
  "\n### 发展阶段分析（基于修改后标准）\n",
  "本研究将认知发展分为四个关键阶段，现使用修改后的显著性标准重新分析：\n",
  "1. **早期快速发展阶段 (6-8岁)**: 基础执行功能快速建立期\n",
  "2. **中期持续发展阶段 (9-12岁)**: 执行功能稳定提升期\n",
  "3. **后期发展阶段 (13-15岁)**: 青春期执行功能调整期\n",
  "4. **稳定阶段 (16-18岁)**: 执行功能成熟期\n\n"
)

# 添加发展完成度分析总结
report_content <- paste0(
  report_content,
  "### 发展完成度分析（基于修改后标准）\n",
  "使用修改后的显著性标准，重新分析了在12、13、14、15、16岁等关键年龄节点的发展完成程度。\n\n",
  
  "### d_value指标的重要性\n",
  "d'值作为信号检测理论的核心指标，具有以下优势：\n",
  "- **敏感性更高**: 比准确率更能反映细微的抑制控制能力变化\n",
  "- **理论基础**: 基于信号检测理论，区分辨别力和反应偏好\n",
  "- **发展敏感**: 能够捕捉儿童青少年期抑制控制的发展轨迹\n\n",
  
  "### 方法学改进与意义\n",
  "本次修改的主要特点：\n",
  "- **更严格的显著性标准**: 要求 |derivatives_6mo_increment| ≥ ", significance_threshold, "\n",
  "- **增加核心指标**: d'值提供更精确的抑制控制测量\n",
  "- **减少假阳性**: 通过绝对值阈值降低微小变化被错误识别为显著的风险\n",
  "- **更保守的发展期识别**: 只有达到实质性变化幅度的时期才被认为是显著发展期\n",
  "- **提高临床意义**: 识别的显著发展期更可能具有实际的生理或心理学意义\n\n"
)

# 写入报告文件
writeLines(report_content, file.path(OUTPUT_PATH, "modified_gam_analysis_report_with_dvalue.md"))
cat("生成修改后分析报告（包含d_value）:", file.path(OUTPUT_PATH, "modified_gam_analysis_report_with_dvalue.md"), "\n")

# 打印修改后平滑项显著性汇总表
cat("\n=== 修改后平滑项显著性汇总表 ===\n")
cat(rep("=", 90), "\n", sep = "")
print(significance_summary[, c("Metric", "F_Value", "P_Value", "Significance_Level",
                               "EDF", "R_squared", "Modified_Sig_Count", "Modified_Sig_Pct")], row.names = FALSE)

# 打印方法比较表
cat("\n=== 方法比较汇总 ===\n")
cat("显著性标准: |derivatives_6mo_increment| ≥", significance_threshold, "\n")
cat(rep("=", 80), "\n", sep = "")
print(method_comparison, row.names = FALSE)

# 特别输出核心抑制指标结果
core_inhibition_found <- successful_metrics[grepl("Stroop_Incongruent|GoNoGo_Nogo_ACC|GoNoGo_d_value", successful_metrics)]
if (length(core_inhibition_found) > 0) {
  cat("\n=== 核心抑制指标结果总结 ===\n")
  cat("检测到的核心抑制指标:", paste(core_inhibition_found, collapse = ", "), "\n")
  cat(rep("=", 80), "\n", sep = "")
  
  core_summary <- summary_stats_display[summary_stats_display$Metric %in% core_inhibition_found, 
                                        c("Metric", "Task", "Variable", "AIC", "R_squared", 
                                          "Smooth_F_Value", "Smooth_P_Value", "Significant_Periods_Count", 
                                          "Significant_Periods_Pct")]
  print(core_summary, row.names = FALSE)
}

# 清理并行计算资源
stopCluster(cl)
cat("\n=== 修改后GAM分析完成（包含d_value） ===\n")
cat("所有结果已保存到", OUTPUT_PATH, "\n")
cat("使用的显著性阈值:", significance_threshold, "\n")
cat("分析的指标数量:", length(successful_metrics), "\n")
if (length(core_inhibition_found) > 0) {
  cat("核心抑制指标:", paste(core_inhibition_found, collapse = ", "), "\n")
}

# 显示最终汇总
cat("\n最终修改后GAM结果汇总（包含d_value）:\n")
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(summary_stats_display[, c("Metric", "AIC", "R_squared", "Smooth_F_Value", 
                                               "Smooth_P_Value", "Smooth_EDF", "Significant_Periods_Count",
                                               "Significant_Periods_Pct", "Significance_Threshold")], 
                     digits = 3))
} else {
  print(summary_stats_display[, c("Metric", "AIC", "R_squared", "Smooth_F_Value", 
                                  "Smooth_P_Value", "Smooth_EDF", "Significant_Periods_Count",
                                  "Significant_Periods_Pct", "Significance_Threshold")])
}