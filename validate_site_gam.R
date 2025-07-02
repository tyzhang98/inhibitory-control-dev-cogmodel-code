# 加载必要的包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  mgcv,           # GAM建模
  readxl,         # 读取Excel文件
  dplyr,          # 数据处理
  tidyr,          # 数据整理
  ggplot2,        # 绘图
  viridis,        # 颜色
  patchwork,      # 图形组合
  openxlsx,       # 写Excel文件
  gridExtra,      # 图形布局
  RColorBrewer,   # 颜色方案
  scales,         # 标度
  broom,          # 模型结果整理
  car,            # 方差分析
  plotly,         # 交互图形
  knitr,          # 表格
  VIM,            # 缺失值处理
  outliers,       # 异常值检测
  boot,           # Bootstrap
  moments         # 偏度峰度计算
)

# 设置图形参数
theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  axis.title = element_text(size = 14, face = "bold"),
                  strip.text = element_text(size = 12, face = "bold"),
                  legend.title = element_text(size = 12, face = "bold")))

# 数据加载和预处理函数
load_and_prepare_data <- function(stroop_path = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Stroop-baseline.xlsx",
                                  gonogo_path = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Gonogo-baseline.xlsx") {
  
  cat("=== 加载数据文件 ===\n")
  # 加载数据
  tryCatch({
    df_stroop <- read_excel(stroop_path)
    df_gonogo <- read_excel(gonogo_path)
    cat("数据文件加载成功\n")
  }, error = function(e) {
    stop(paste("加载文件错误:", e$message))
  })
  
  # 处理Stroop数据
  stroop_metrics <- c('Incongruent', 'Congruent', 'Neutral', 'Interference')
  available_stroop <- intersect(stroop_metrics, colnames(df_stroop))
  
  cat("使用的Stroop指标:", paste(available_stroop, collapse = ", "), "\n")
  
  # 处理GoNoGo数据（添加新的核心指标）
  gonogo_metrics <- c('Go_ACC', 'Nogo_ACC', 'Go_RT_ms', "Gonogo_d'")  # 添加了Gonogo_d'
  available_gonogo <- intersect(gonogo_metrics, colnames(df_gonogo))
  
  cat("使用的GoNoGo指标:", paste(available_gonogo, collapse = ", "), "\n")
  
  # 创建综合数据集
  cat("\n=== 创建综合数据集 ===\n")
  all_data <- list()
  
  # 处理Stroop数据
  for (metric in available_stroop) {
    if (metric %in% colnames(df_stroop)) {
      temp_df <- df_stroop %>%
        select(Subject, Age, all_of(metric)) %>%
        rename(value = all_of(metric)) %>%
        mutate(
          task = "Stroop",
          metric = metric,
          task_metric = paste0("Stroop_", metric)
        )
      all_data[[paste0("Stroop_", metric)]] <- temp_df
    }
  }
  
  # 处理GoNoGo数据
  for (metric in available_gonogo) {
    if (metric %in% colnames(df_gonogo)) {
      temp_df <- df_gonogo %>%
        select(Subject, Age, all_of(metric)) %>%
        rename(value = all_of(metric)) %>%
        mutate(
          task = "GoNoGo",
          metric = metric,
          task_metric = paste0("GoNoGo_", metric)
        )
      all_data[[paste0("GoNoGo_", metric)]] <- temp_df
    }
  }
  
  # 合并所有数据
  df_combined <- bind_rows(all_data) %>%
    rename(age = Age) %>%
    filter(!is.na(value), !is.na(age))
  
  cat("合并数据集:", nrow(df_combined), "个观测值，覆盖",
      length(unique(df_combined$task_metric)), "个指标\n")
  
  # Z分数标准化
  df_combined <- df_combined %>%
    group_by(task_metric) %>%
    mutate(value_z = scale(value)[,1]) %>%
    ungroup()
  
  return(df_combined)
}

# 异常值检测和清理函数
clean_outliers <- function(data, z_threshold = 4) {
  cat("原始样本量:", nrow(data), "\n")
  
  # 单变量异常值检测
  data_clean <- data %>%
    group_by(task_metric) %>%
    filter(abs(value_z) <= z_threshold) %>%
    ungroup()
  
  cat("单变量异常值清理后:", nrow(data_clean), "\n")
  
  # 多变量异常值检测（马氏距离）
  if (nrow(data_clean) > 100) {
    # 计算每个任务指标的马氏距离
    outlier_subjects <- c()
    
    for (tm in unique(data_clean$task_metric)) {
      subset_data <- data_clean %>% filter(task_metric == tm)
      if (nrow(subset_data) > 10) {
        # 计算马氏距离
        numeric_vars <- subset_data %>% select(age, value_z)
        if (nrow(numeric_vars) > ncol(numeric_vars)) {
          tryCatch({
            mahal_dist <- mahalanobis(numeric_vars,
                                      colMeans(numeric_vars),
                                      cov(numeric_vars))
            threshold <- qchisq(0.975, ncol(numeric_vars))
            outliers <- subset_data$Subject[mahal_dist > threshold]
            outlier_subjects <- c(outlier_subjects, outliers)
          }, error = function(e) {
            cat("马氏距离计算失败，任务:", tm, "\n")
          })
        }
      }
    }
    
    # 移除异常值
    data_final <- data_clean %>%
      filter(!Subject %in% unique(outlier_subjects))
    
    cat("多变量异常值清理后:", nrow(data_final), "\n")
  } else {
    data_final <- data_clean
  }
  
  return(data_final)
}

# GAM模型拟合函数
fit_gam_model <- function(data, outcome_var = "value_z", age_var = "age") {
  
  # 构建GAM模型公式
  formula_str <- paste(outcome_var, "~ s(", age_var, ", bs='cr')")
  gam_formula <- as.formula(formula_str)
  
  # 拟合GAM模型
  tryCatch({
    gam_model <- gam(gam_formula, data = data, method = "REML")
    return(gam_model)
  }, error = function(e) {
    cat("GAM拟合失败:", e$message, "\n")
    return(NULL)
  })
}

# 改进的GAM导数计算函数
compute_gam_derivatives <- function(gam_model, data, age_range = NULL,
                                    n_bootstrap = 3000) {  # 新增：最小效应量阈值
  if (is.null(gam_model)) {
    return(NULL)
  }
  
  if (is.null(age_range)) {
    age_range <- seq(min(data$age, na.rm = TRUE),
                     max(data$age, na.rm = TRUE), by = 0.1)
  }
  
  # 创建预测数据
  pred_data <- data.frame(age = age_range)
  
  # 获取预测值和标准误差
  predictions <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
  fitted_values <- predictions$fit
  se_fit <- predictions$se.fit
  
  # 计算95%置信区间
  pred_lower <- fitted_values - 1.96 * se_fit
  pred_upper <- fitted_values + 1.96 * se_fit
  
  # 计算一阶导数
  # 使用finite differences方法
  h <- age_range[2] - age_range[1]
  derivatives <- rep(NA, length(fitted_values))
  
  # 中心差分法计算导数
  for (i in 2:(length(fitted_values) - 1)) {
    derivatives[i] <- (fitted_values[i + 1] - fitted_values[i - 1]) / (2 * h)
  }
  
  # 边界点使用前向/后向差分
  derivatives[1] <- (fitted_values[2] - fitted_values[1]) / h
  derivatives[length(derivatives)] <- 
    (fitted_values[length(fitted_values)] - fitted_values[length(fitted_values) - 1]) / h
  
  # Bootstrap置信区间for导数
  cat("    计算导数的bootstrap置信区间 (", n_bootstrap, "次抽样)...\n")
  
  n_obs <- nrow(data)
  n_bootstrap_actual <- min(n_bootstrap, 3000)
  
  if (n_obs < 30) {
    n_bootstrap_actual <- min(n_bootstrap_actual, 500)
  }
  
  deriv_bootstrap <- matrix(NA, nrow = n_bootstrap_actual, ncol = length(age_range))
  
  for (i in 1:n_bootstrap_actual) {
    if (i %% 200 == 0) {
      cat("    Bootstrap进度:", i, "/", n_bootstrap_actual, "\n")
    }
    
    # Bootstrap抽样
    boot_indices <- sample(1:n_obs, size = n_obs, replace = TRUE)
    boot_data <- data[boot_indices, ]
    
    tryCatch({
      # 拟合bootstrap模型
      boot_gam <- gam(value_z ~ s(age, bs = 'cr'),
                      data = boot_data, method = "REML")
      
      # 预测
      boot_pred <- predict(boot_gam, newdata = pred_data)
      
      # 计算导数
      boot_deriv <- rep(NA, length(boot_pred))
      for (j in 2:(length(boot_pred) - 1)) {
        boot_deriv[j] <- (boot_pred[j + 1] - boot_pred[j - 1]) / (2 * h)
      }
      boot_deriv[1] <- (boot_pred[2] - boot_pred[1]) / h
      boot_deriv[length(boot_deriv)] <- 
        (boot_pred[length(boot_pred)] - boot_pred[length(boot_pred) - 1]) / h
      
      deriv_bootstrap[i, ] <- boot_deriv
      
    }, error = function(e) {
      # 如果bootstrap失败，跳过这次迭代
    })
  }
  
  # 计算导数置信区间
  valid_bootstrap <- !is.na(deriv_bootstrap[, 1])
  if (sum(valid_bootstrap) > 10) {
    deriv_lower <- apply(deriv_bootstrap[valid_bootstrap, ], 2,
                         quantile, probs = 0.025, na.rm = TRUE)
    deriv_upper <- apply(deriv_bootstrap[valid_bootstrap, ], 2,
                         quantile, probs = 0.975, na.rm = TRUE)
  } else {
    # 如果bootstrap失败，使用简单的标准误差估计
    deriv_se <- sd(derivatives, na.rm = TRUE)
    deriv_lower <- derivatives - 1.96 * deriv_se
    deriv_upper <- derivatives + 1.96 * deriv_se
  }
  
  # 显著性检验（加入效应量阈值判断）
  significant <- (deriv_lower > 0) | (deriv_upper < 0)
  
  # 新增：效应量阈值筛选 - 只有绝对值大于最小效应量阈值的才认为是显著的
  effect_size_mask <- abs(derivatives) >= min_effect_size
  significant <- significant & effect_size_mask
  
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
  
  # 模型统计信息
  model_summary <- summary(gam_model)
  
  # 提取平滑项显著性检验结果
  smooth_test <- model_summary$s.table
  smooth_f_value <- smooth_test[1, "F"]  # F统计量
  smooth_p_value <- smooth_test[1, "p-value"]  # p值
  smooth_ref_df <- smooth_test[1, "Ref.df"]  # 参考自由度
  
  results <- list(
    age = age_range,
    predictions = fitted_values,
    pred_lower = pred_lower,
    pred_upper = pred_upper,
    derivatives = derivatives,
    deriv_lower = deriv_lower,
    deriv_upper = deriv_upper,
    significant = significant_cleaned,
    aic = AIC(gam_model),
    deviance_explained = model_summary$dev.expl,
    edf = sum(model_summary$edf),
    gcv_score = gam_model$gcv.ubre,
    residual_se = sqrt(gam_model$sig2),
    n_bootstrap_actual = sum(valid_bootstrap),
    smooth_f_value = smooth_f_value,       # 新增：F统计量
    smooth_p_value = smooth_p_value,       # 新增：p值
    smooth_ref_df = smooth_ref_df,         # 新增：参考自由度
    model = gam_model,
    model_summary = model_summary
  )
  
  return(results)
}

# 倒U型模式检测函数
detect_inverted_u_pattern <- function(results, age_range,
                                      min_peak_age = 6, max_peak_age = 15) {
  
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

# 创建可视化函数（修改版 - 添加外边框和显著发展年龄标记）
create_comprehensive_plots <- function(gam_results, successful_metrics, data_final) {
  cat("\n=== 创建可视化图形 ===\n")
  
  # 1. 所有轨迹图（修改部分）
  plot_data <- data.frame()
  
  for (task_metric in successful_metrics) {
    results <- gam_results[[task_metric]]
    if (!is.null(results)) {
      temp_data <- data.frame(
        age = results$age,
        predictions = results$predictions,
        pred_lower = results$pred_lower,
        pred_upper = results$pred_upper,
        significant = results$significant,  # 添加显著性信息
        task_metric = task_metric
      )
      plot_data <- rbind(plot_data, temp_data)
    }
  }
  
  # 颜色和线型设置
  plot_data$task <- ifelse(grepl("Stroop", plot_data$task_metric), "Stroop", "GoNoGo")
  
  # 替换task_metric中的下划线为点号，用于图例显示
  plot_data$task_metric_display <- gsub("_", " ", plot_data$task_metric)
  
  # 计算Stroop不一致和NoGo准确率的最晚显著发展年龄
  calculate_specific_latest_age <- function(metric_name) {
    if (metric_name %in% successful_metrics) {
      results <- gam_results[[metric_name]]
      if (!is.null(results)) {
        sig_ages <- results$age[results$significant]
        if (length(sig_ages) > 0) {
          return(max(sig_ages))
        }
      }
    }
    return(NA)
  }
  
  # 获取特定指标的最晚显著发展年龄
  stroop_incongruent_latest <- calculate_specific_latest_age("Stroop_Incongruent")
  nogo_acc_latest <- calculate_specific_latest_age("GoNoGo_Nogo_ACC")
  
  cat("Stroop不一致最晚显著发展年龄:", stroop_incongruent_latest, "\n")
  cat("NoGo准确率最晚显著发展年龄:", nogo_acc_latest, "\n")
  
  p1 <- ggplot(plot_data, aes(x = age, y = predictions)) +
    # 置信区间
    geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = task),
                alpha = 0.2) +
    # 主轨迹线
    geom_line(aes(color = task, linetype = task_metric_display), linewidth = 1.2) +
    # 添加显著发展期的红色圆点
    geom_point(data = plot_data[plot_data$significant, ],
               aes(x = age, y = predictions),
               color = "red", size = 0.6, alpha = 0.9, shape = 16) +
    
    # 添加Stroop不一致和NoGo准确率的最晚显著发展年龄竖线和填充区域
    {if (!is.na(stroop_incongruent_latest) && !is.na(nogo_acc_latest)) {
      # 中间区域填充
      list(
        annotate("rect",
                 xmin = min(stroop_incongruent_latest, nogo_acc_latest, na.rm = TRUE),
                 xmax = max(stroop_incongruent_latest, nogo_acc_latest, na.rm = TRUE),
                 ymin = -Inf, ymax = Inf,
                 alpha = 0.29, fill = "#3690c0"),
        # Stroop不一致最晚年龄竖线
        geom_vline(xintercept = stroop_incongruent_latest,
                   color = "black", linetype = "dashed",
                   linewidth = 0.5, alpha = 0.8),
        # NoGo准确率最晚年龄竖线  
        geom_vline(xintercept = nogo_acc_latest,
                   color = "black",  linetype = "dashed",
                   linewidth = 0.5, alpha = 0.8),
        # 添加文本标签
        annotate("text", x = stroop_incongruent_latest, y = 1.8,
                 label = paste0(""),
                 color = "#0868ac", fontface = "bold", size = 3.5,
                 hjust = ifelse(stroop_incongruent_latest > nogo_acc_latest, 1.1, -0.1)),
        annotate("text", x = nogo_acc_latest, y = 1.6,
                 label = paste0(""),
                 color = "#006d2c", fontface = "bold", size = 3.5,
                 hjust = ifelse(nogo_acc_latest > stroop_incongruent_latest, 1.1, -0.1))
      )
    } else if (!is.na(stroop_incongruent_latest)) {
      # 仅显示Stroop不一致竖线
      list(
        geom_vline(xintercept = stroop_incongruent_latest,
                   color = "black", linetype = "dashed",
                   linewidth = 0.5, alpha = 0.8),
        annotate("text", x = stroop_incongruent_latest, y = 1.8,
                 
                 color = "#0868ac", fontface = "bold", size = 3.5,
                 hjust = 0.5)
      )
    } else if (!is.na(nogo_acc_latest)) {
      # 仅显示NoGo准确率竖线
      list(
        geom_vline(xintercept = nogo_acc_latest,
                   color = "black", linetype = "dashed",
                   linewidth = 0.5, alpha = 0.8),
        annotate("text", x = nogo_acc_latest, y = 1.6,
                 
                 color = "#006d2c", fontface = "bold", size = 3.5,
                 hjust = 0.5)
      )
    } else {
      list()
    }} +
    
    scale_color_manual(values = c("Stroop" = "#0868ac", "GoNoGo" = "#006d2c")) +
    scale_fill_manual(values = c("Stroop" = "#0868ac", "GoNoGo" = "#006d2c")) +
    labs(x = "Age (years)", y = "Z-scored Performance",
         title = "") +
    
    theme_minimal() +
    theme(
      legend.position = "right",
      # 图例字体调小
      legend.text = element_text(size = 6),
      # 删除所有图例标题
      legend.title = element_blank(),
      # 调整图例间距和尺寸
      legend.key.size = unit(0.4, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.35, "cm"),
      legend.spacing.y = unit(0.12, "cm"),
      legend.spacing.x = unit(0.12, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 5, unit = "pt"),
      legend.box.spacing = unit(0.12, "cm"),
      # 添加外边框加粗
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      # 其他样式保持不变
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 12, face = "bold")
    ) +
    xlim(6, 19) +
    ylim(-2, 2) +
    scale_x_continuous(breaks = seq(6, 19, 1))
  
  # 2. 核心指标单独图形（更新了核心指标列表）
  core_metrics <- c("Stroop_Incongruent", "Stroop_Interference", "GoNoGo_Nogo_ACC", "GoNoGo_Gonogo_d'")
  available_core <- intersect(core_metrics, successful_metrics)
  
  # 调试信息：检查核心指标匹配情况
  cat("期望的核心指标:", paste(core_metrics, collapse = ", "), "\n")
  cat("成功分析的指标:", paste(successful_metrics, collapse = ", "), "\n")
  cat("匹配的核心指标:", paste(available_core, collapse = ", "), "\n")
  
  core_plots <- list()
  
  if (length(available_core) > 0) {
    cat("\n创建核心指标单独图形...\n")
    
    for (core_metric in available_core) {
      cat("  创建", core_metric, "图形\n")
      
      # 获取GAM拟合结果
      results <- gam_results[[core_metric]]
      if (is.null(results)) next
      
      # 获取原始数据点
      raw_data <- data_final[data_final$task_metric == core_metric, ]
      
      # 创建拟合曲线数据
      fit_data <- data.frame(
        age = results$age,
        predictions = results$predictions,
        pred_lower = results$pred_lower,
        pred_upper = results$pred_upper,
        significant = results$significant  # 添加显著性信息
      )
      
      # 为原始数据计算对应的拟合值（用于连线）
      raw_data$fitted_values <- predict(results$model, newdata = raw_data)
      
      # 确定任务类型和颜色
      task_type <- ifelse(grepl("Stroop", core_metric), "Stroop", "GoNoGo")
      main_color <- ifelse(task_type == "Stroop", "#0868ac", "#006d2c")  # 蓝色 vs 绿色
      fitted_color <- ifelse(task_type == "Stroop", "#4292c6", "#41ab5d")  # 浅蓝色 vs 浅绿色
      
      # 为颜色映射选择合适的颜色方案
      color_scheme <- ifelse(task_type == "Stroop", "Blues", "Greens")
      
      # 创建图形
      p_core <- ggplot() +
        # 置信区间（最底层）
        geom_ribbon(data = fit_data, aes(x = age, ymin = pred_lower, ymax = pred_upper),
                    alpha = 0.2, fill = main_color) +
        # 原始数据点到拟合值的连线（残差线）
        geom_segment(data = raw_data,
                     aes(x = age, y = value_z, xend = age, yend = fitted_values),
                     alpha = 0.5, color = "grey60", linewidth = 0.2) +
        # 拟合曲线（主曲线）
        geom_line(data = fit_data, aes(x = age, y = predictions),
                  color = main_color, linewidth = 2.5, alpha = 0.8) +
        # 拟合值散点（对应原始数据点的拟合位置）
        geom_point(data = raw_data, aes(x = age, y = fitted_values),
                   color = fitted_color, size = 2, alpha = 0.8, shape = 17) +  # 三角形
        # 原始数据点（颜色映射到Z值）
        geom_point(data = raw_data, aes(x = age, y = value_z, color = value_z),
                   size = 2.5, alpha = 0.9, shape = 16) +  # 圆形
        # Z值颜色映射
        scale_color_distiller(type = "seq", palette = color_scheme, direction = 1,
                              name = "Z-score",
                              guide = guide_colorbar(title.position = "top",
                                                     barwidth = 8, barheight = 0.8)) +
        # 针状图 - X轴上的原始数据分布
        geom_rug(data = raw_data, aes(x = age),
                 color = main_color, alpha = 0.6, length = unit(0.02, "npc")) +
        # 针状图 - Y轴上的原始数据分布
        geom_rug(data = raw_data, aes(y = value_z),
                 color = main_color, alpha = 0.6, length = unit(0.02, "npc")) +
        # 添加显著发展期标记 - 改为红色圆点
        geom_point(data = fit_data[fit_data$significant, ],
                   aes(x = age, y = predictions),
                   color = "red", size = 1.5, alpha = 0.9, shape = 16) +  # 红色圆点
        # 格式设置
        labs(x = "Age (years)",
             y = "Z-scored Performance",
             subtitle = paste0("",
                               "")) +
        scale_x_continuous(breaks = seq(6, 19, 1), limits = c(6, 19)) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey50"),
          # X、Y轴标题和刻度字体调大
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),

          panel.grid = element_blank(),      # 隐藏所有网格线
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),

          legend.position = "bottom",  # 将图例位置改为顶部
          legend.direction = "horizontal",
          legend.margin = margin(b = 0.1)  # 在图例下方添加一些边距
        )
      
      core_plots[[core_metric]] <- p_core
    }
  }
  
  return(list(
    trajectory_plot = p1,
    core_plots = core_plots
  ))
}




# 结果打印函数
print_comprehensive_results <- function(gam_results, inverted_u_stats,
                                        stage_stats, completion_stats,
                                        successful_metrics) {
  
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
  
  # 发展阶段变化率分析
  cat("\n=== 发展阶段变化率比较分析（核心抑制控制指标）===\n")
  cat(rep("=", 60), "\n", sep = "")
  
  core_inhibition_metrics <- successful_metrics[
    grepl("Stroop_Incongruent|GoNoGo_Nogo_ACC|GoNoGo_Gonogo_d'|Gonogo_d", successful_metrics)
  ]
  
  if (length(core_inhibition_metrics) > 0) {
    cat("分析的核心抑制控制指标: ", paste(core_inhibition_metrics, collapse = ", "), "\n")
    
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
    cat("未找到核心抑制控制指标\n")
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

# 保存结果函数
save_comprehensive_results <- function(gam_results, data_final, successful_metrics,
                                       summary_stats, inverted_u_stats, stage_stats,
                                       completion_stats, plots, output_dir = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/output") {
  
  cat("\n=== 保存分析结果 ===\n")
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 准备详细结果数据框
  age_range <- seq(min(data_final$age), max(data_final$age), by = 0.1)
  detailed_results <- data.frame(age = age_range)
  
  for (task_metric in successful_metrics) {
    results <- gam_results[[task_metric]]
    if (!is.null(results)) {
      # 插值到统一的年龄网格
      detailed_results[[paste0(task_metric, "_prediction")]] <-
        approx(results$age, results$predictions, xout = age_range)$y
      detailed_results[[paste0(task_metric, "_derivative")]] <-
        approx(results$age, results$derivatives, xout = age_range)$y
      detailed_results[[paste0(task_metric, "_significant")]] <-
        approx(results$age, as.numeric(results$significant), xout = age_range)$y
      detailed_results[[paste0(task_metric, "_deriv_lower")]] <-
        approx(results$age, results$deriv_lower, xout = age_range)$y
      detailed_results[[paste0(task_metric, "_deriv_upper")]] <-
        approx(results$age, results$deriv_upper, xout = age_range)$y
    }
  }
  
  # 准备热力图数据（更细的网格，用于热力图绘制）
  age_min <- min(data_final$age, na.rm = TRUE)
  age_max <- max(data_final$age, na.rm = TRUE)
  heatmap_age_seq <- seq(age_min, age_max, by = 0.05)  # 更细的网格
  heatmap_data_export <- data.frame(age = heatmap_age_seq)
  
  cat("准备热力图数据导出...\n")
  
  for (task_metric in successful_metrics) {
    results <- gam_results[[task_metric]]
    if (!is.null(results)) {
      cat("  处理", task_metric, "热力图数据\n")
      
      # 插值到热力图网格
      pred_interp <- approx(results$age, results$predictions, xout = heatmap_age_seq, rule = 2)$y
      deriv_interp <- approx(results$age, results$derivatives, xout = heatmap_age_seq, rule = 2)$y
      sig_interp <- approx(results$age, as.numeric(results$significant), xout = heatmap_age_seq, rule = 2)$y
      deriv_lower_interp <- approx(results$age, results$deriv_lower, xout = heatmap_age_seq, rule = 2)$y
      deriv_upper_interp <- approx(results$age, results$deriv_upper, xout = heatmap_age_seq, rule = 2)$y
      
      # 添加到热力图数据框
      heatmap_data_export[[paste0(task_metric, "_prediction")]] <- pred_interp
      heatmap_data_export[[paste0(task_metric, "_derivative")]] <- deriv_interp
      heatmap_data_export[[paste0(task_metric, "_significant")]] <- sig_interp >= 0.5  # 转换为逻辑值
      heatmap_data_export[[paste0(task_metric, "_deriv_lower")]] <- deriv_lower_interp
      heatmap_data_export[[paste0(task_metric, "_deriv_upper")]] <- deriv_upper_interp
    }
  }
  
  # 准备汇总统计
  summary_df <- data.frame()
  
  for (task_metric in successful_metrics) {
    results <- gam_results[[task_metric]]
    if (!is.null(results)) {
      
      metric_data <- data_final[data_final$task_metric == task_metric, ]
      
      # 基本统计
      stats_row <- data.frame(
        Metric = task_metric,
        Task = ifelse(grepl("Stroop", task_metric), "Stroop", "GoNoGo"),
        Variable = gsub("^(Stroop_|GoNoGo_)", "", task_metric),
        N_Observations = nrow(metric_data),
        AIC = results$aic,
        R_squared = results$deviance_explained,
        EDF = results$edf,
        GCV_Score = results$gcv_score,
        Residual_SE = results$residual_se,
        Bootstrap_Samples = results$n_bootstrap_actual,
        Smooth_F_Value = results$smooth_f_value,      # 新增：F统计量
        Smooth_P_Value = results$smooth_p_value,      # 新增：p值
        Smooth_Ref_DF = results$smooth_ref_df         # 新增：参考自由度
      )
      
      # 显著性统计
      sig_periods <- results$significant
      stats_row$Significant_Periods_Pct <- sum(sig_periods) / length(sig_periods) * 100
      stats_row$Max_Derivative <- max(abs(results$derivatives), na.rm = TRUE)
      stats_row$Mean_Derivative <- mean(results$derivatives, na.rm = TRUE)
      
      if (any(sig_periods)) {
        sig_ages <- results$age[sig_periods]
        stats_row$Sig_Age_Start <- min(sig_ages)
        stats_row$Sig_Age_End <- max(sig_ages)
      } else {
        stats_row$Sig_Age_Start <- NA
        stats_row$Sig_Age_End <- NA
      }
      
      max_change_idx <- which.max(abs(results$derivatives))
      stats_row$Peak_Change_Age <- results$age[max_change_idx]
      
      # 倒U型统计
      u_stats <- inverted_u_stats[[task_metric]]
      stats_row$Has_Inverted_U <- u_stats$has_inverted_u
      stats_row$U_Peak_Age <- u_stats$peak_age
      stats_row$U_Peak_Value <- u_stats$peak_value
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
      stats_row$Sig_Duration <- completion_info$sig_duration
      cutoff_ages <- c(12, 13, 14, 15, 16)
      for (cutoff in cutoff_ages) {
        completion_key <- paste0("completion_", cutoff)
        duration_key <- paste0("pre_", cutoff, "_duration")
        if (completion_key %in% names(completion_info)) {
          stats_row[[paste0("Completion_", cutoff, "_Pct")]] <- completion_info[[completion_key]]
          stats_row[[paste0("Pre_", cutoff, "_Duration")]] <- completion_info[[duration_key]]
        }
      }
      
      summary_df <- rbind(summary_df, stats_row)
    }
  }
  
  # 转换列表为数据框
  inverted_u_df <- do.call(rbind, lapply(names(inverted_u_stats), function(metric) {
    stats <- inverted_u_stats[[metric]]
    data.frame(
      Metric = metric,
      Has_Inverted_U = stats$has_inverted_u,
      Peak_Age = ifelse(is.na(stats$peak_age), NA, stats$peak_age),
      Peak_Value = ifelse(is.na(stats$peak_value), NA, stats$peak_value),
      Increase_Duration = stats$increase_duration,
      Decrease_Duration = stats$decrease_duration,
      stringsAsFactors = FALSE
    )
  }))
  
  stage_stats_df <- do.call(rbind, lapply(names(stage_stats), function(metric) {
    stats <- stage_stats[[metric]]
    data.frame(
      Metric = metric,
      Early_6_8_Mean_Rate = stats$early_rapid_6_8_mean_rate,
      Early_6_8_Abs_Rate = stats$early_rapid_6_8_abs_mean_rate,
      Early_6_8_Duration = stats$early_rapid_6_8_duration,
      Early_6_8_Age_Range = stats$early_rapid_6_8_age_range,
      Middle_9_12_Mean_Rate = stats$middle_sustained_9_12_mean_rate,
      Middle_9_12_Abs_Rate = stats$middle_sustained_9_12_abs_mean_rate,
      Middle_9_12_Duration = stats$middle_sustained_9_12_duration,
      Middle_9_12_Age_Range = stats$middle_sustained_9_12_age_range,
      Late_13_15_Mean_Rate = stats$late_plateau_13_15_mean_rate,
      Late_13_15_Abs_Rate = stats$late_plateau_13_15_abs_mean_rate,
      Late_13_15_Duration = stats$late_plateau_13_15_duration,
      Late_13_15_Age_Range = stats$late_plateau_13_15_age_range,
      Final_16_18_Mean_Rate = stats$final_stable_16_18_mean_rate,
      Final_16_18_Abs_Rate = stats$final_stable_16_18_abs_mean_rate,
      Final_16_18_Duration = stats$final_stable_16_18_duration,
      Final_16_18_Age_Range = stats$final_stable_16_18_age_range,
      stringsAsFactors = FALSE
    )
  }))
  
  completion_df <- do.call(rbind, lapply(names(completion_stats), function(metric) {
    stats <- completion_stats[[metric]]
    row_data <- data.frame(Metric = metric, stringsAsFactors = FALSE)
    for (name in names(stats)) {
      row_data[[name]] <- stats[[name]]
    }
    return(row_data)
  }))
  
  # 保存到Excel文件
  tryCatch({
    output_path <- file.path(output_dir, "comprehensive_gam_analysis_results.xlsx")
    
    # 创建工作簿
    wb <- createWorkbook()
    
    # 添加工作表
    addWorksheet(wb, "Clean_Data")
    addWorksheet(wb, "GAM_Results")
    addWorksheet(wb, "Heatmap_Data")  # 新增：热力图数据工作表
    addWorksheet(wb, "Summary_Statistics")
    addWorksheet(wb, "Inverted_U_Analysis")
    addWorksheet(wb, "Development_Stage_Rates")
    addWorksheet(wb, "Development_Completion")
    
    # 写入数据
    writeData(wb, "Clean_Data", data_final)
    writeData(wb, "GAM_Results", detailed_results)
    writeData(wb, "Heatmap_Data", heatmap_data_export)  # 新增：写入热力图数据
    writeData(wb, "Summary_Statistics", summary_df)
    writeData(wb, "Inverted_U_Analysis", inverted_u_df)
    writeData(wb, "Development_Stage_Rates", stage_stats_df)
    writeData(wb, "Development_Completion", completion_df)
    
    # 保存文件
    saveWorkbook(wb, output_path, overwrite = TRUE)
    cat("结果已保存到:", output_path, "\n")
    
  }, error = function(e) {
    cat("保存Excel文件失败:", e$message, "\n")
    # 备用：保存为CSV文件
    write.csv(summary_df, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
    write.csv(detailed_results, file.path(output_dir, "gam_results.csv"), row.names = FALSE)
    write.csv(heatmap_data_export, file.path(output_dir, "heatmap_data.csv"), row.names = FALSE)  # 新增：保存热力图数据为CSV
    cat("结果已保存为CSV格式\n")
  })
  
  # 保存图形
  tryCatch({
    ggsave(file.path(output_dir, "all_metrics_gam_fits.png"),
           plots$trajectory_plot, width = 5.5, height = 3.4, dpi = 800)
    
    # 保存核心指标单独图形
    if (!is.null(plots$core_plots) && length(plots$core_plots) > 0) {
      core_dir <- file.path(output_dir, "core_metrics_plots")
      if (!dir.exists(core_dir)) {
        dir.create(core_dir, recursive = TRUE)
      }
      
      for (metric_name in names(plots$core_plots)) {
        filename <- file.path(core_dir, paste0(metric_name, "_individual_fit.png"))
        ggsave(filename, plots$core_plots[[metric_name]],
               width = 2, height = 1.5, dpi = 800)
        cat("核心指标图形已保存:", filename, "\n")
      }
    }
    
    cat("图形已保存\n")
  }, error = function(e) {
    cat("保存图形失败:", e$message, "\n")
  })
  
  return(summary_df)
}

# 主分析函数
main_gam_analysis <- function(stroop_path = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Stroop-baseline.xlsx",
                              gonogo_path = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Gonogo-baseline.xlsx") {
  
  cat("开始GAM分析流程...\n")
  
  # 1. 数据加载和预处理
  data_combined <- load_and_prepare_data(stroop_path, gonogo_path)
  
  # 2. 数据清理
  data_final <- clean_outliers(data_combined)
  
  # 3. 定义年龄范围
  age_range <- seq(min(data_final$age), max(data_final$age), by = 0.1)
  
  # 4. GAM分析
  cat("\n=== GAM分析开始 ===\n")
  gam_results <- list()
  successful_metrics <- c()
  
  for (task_metric in unique(data_final$task_metric)) {
    cat("\n分析", task_metric, "...\n")
    metric_data <- data_final[data_final$task_metric == task_metric, ]
    
    if (nrow(metric_data) < 20) {
      cat("  跳过", task_metric, ": 数据不足 (", nrow(metric_data), "个观测值)\n")
      next
    }
    
    # 拟合GAM模型
    gam_model <- fit_gam_model(metric_data, "value_z", "age")
    
    if (!is.null(gam_model)) {
      # 计算导数和统计
      results <- compute_gam_derivatives(gam_model, metric_data, age_range, n_bootstrap = 3000)
      
      if (!is.null(results)) {
        gam_results[[task_metric]] <- results
        successful_metrics <- c(successful_metrics, task_metric)
        
        cat("  AIC:", round(results$aic, 2),
            ", R²:", round(results$deviance_explained * 100, 1), "%",
            ", EDF:", round(results$edf, 2),
            ", F值:", round(results$smooth_f_value, 2),
            ", p值:", round(results$smooth_p_value, 4), "\n")
      }
    } else {
      cat("  GAM拟合失败\n")
    }
  }
  
  if (length(successful_metrics) == 0) {
    stop("没有成功的GAM分析结果")
  }
  
  cat("\n成功分析:", length(successful_metrics), "个指标\n")
  
  # 5. 倒U型模式检测
  cat("\n=== 检测倒U型发展模式 ===\n")
  inverted_u_stats <- list()
  
  for (task_metric in successful_metrics) {
    results <- gam_results[[task_metric]]
    inverted_u_stats[[task_metric]] <- detect_inverted_u_pattern(results, age_range)
  }
  
  # 6. 发展阶段变化率计算
  cat("\n=== 计算发展阶段变化率 ===\n")
  stage_stats <- list()
  
  for (task_metric in successful_metrics) {
    results <- gam_results[[task_metric]]
    stage_stats[[task_metric]] <- calculate_development_stage_rates(results, age_range)
  }
  
  # 7. 发展完成度计算
  cat("\n=== 计算发展完成度 ===\n")
  completion_stats <- list()
  cutoff_ages <- c(12, 13, 14, 15, 16)
  
  for (task_metric in successful_metrics) {
    results <- gam_results[[task_metric]]
    completion_stats[[task_metric]] <- calculate_development_completion(results, age_range, cutoff_ages)
  }
  
  # 8. 创建可视化
  plots <- create_comprehensive_plots(gam_results, successful_metrics, data_final)
  
  # 9. 打印结果
  print_comprehensive_results(gam_results, inverted_u_stats, stage_stats,
                              completion_stats, successful_metrics)
  
  # 10. 保存结果
  summary_stats <- save_comprehensive_results(gam_results, data_final, successful_metrics,
                                              NULL, inverted_u_stats, stage_stats,
                                              completion_stats, plots)
  
  # 11. 显示图形
  print(plots$trajectory_plot)
  
  # 显示核心指标单独图形
  if (!is.null(plots$core_plots) && length(plots$core_plots) > 0) {
    cat("\n=== 显示核心指标单独图形 ===\n")
    for (metric_name in names(plots$core_plots)) {
      cat("显示", metric_name, "图形\n")
      print(plots$core_plots[[metric_name]])
    }
  }
  
  cat("\n=== GAM分析完成 ===\n")
  
  return(list(
    gam_results = gam_results,
    data_final = data_final,
    successful_metrics = successful_metrics,
    inverted_u_stats = inverted_u_stats,
    stage_stats = stage_stats,
    completion_stats = completion_stats,
    plots = plots,
    summary_stats = summary_stats
  ))
}

# 模型诊断函数
run_model_diagnostics <- function(gam_results, successful_metrics) {
  
  cat("\n=== Model Diagnostic Report ===\n")
  cat(rep("=", 60), "\n", sep = "")
  
  for (metric in successful_metrics) {
    results <- gam_results[[metric]]
    model <- results$model
    
    cat("\n", metric, ":\n", sep = "")
    cat("Model formula: value_z ~ s(age, bs='cr')\n")
    cat("Estimation method: REML\n")
    cat("AIC: ", round(results$aic, 2), "\n")
    cat("Deviance explained: ", round(results$deviance_explained * 100, 1), "%\n")
    cat("Effective degrees of freedom: ", round(results$edf, 2), "\n")
    cat("GCV score: ", round(results$gcv_score, 4), "\n")
    cat("Residual standard error: ", round(results$residual_se, 4), "\n")
    cat("Smooth term F-statistic: ", round(results$smooth_f_value, 2), "\n")
    cat("Smooth term p-value: ", round(results$smooth_p_value, 4), "\n")
    cat("Smooth term reference df: ", round(results$smooth_ref_df, 2), "\n")
    
    # Determine smooth term significance
    if (results$smooth_p_value < 0.001) {
      cat("Smooth term significance: *** (p < 0.001)\n")
    } else if (results$smooth_p_value < 0.01) {
      cat("Smooth term significance: ** (p < 0.01)\n")
    } else if (results$smooth_p_value < 0.05) {
      cat("Smooth term significance: * (p < 0.05)\n")
    } else if (results$smooth_p_value < 0.1) {
      cat("Smooth term significance: . (p < 0.1)\n")
    } else {
      cat("Smooth term significance: Not significant (p >= 0.1)\n")
    }
    
    # Model validation
    tryCatch({
      # Get model residuals
      residuals <- residuals(model)
      
      # Normality test
      if (length(residuals) < 5000) {
        shapiro_test <- shapiro.test(residuals)
        cat("Shapiro-Wilk normality test: W =", round(shapiro_test$statistic, 4),
            ", p =", round(shapiro_test$p.value, 4), "\n")
      }
      
      # Residual statistics
      cat("Residual skewness: ", round(moments::skewness(residuals), 4), "\n")
      cat("Residual kurtosis: ", round(moments::kurtosis(residuals), 4), "\n")
      
    }, error = function(e) {
      cat("Model diagnostic calculation failed\n")
    })
  }
}

# 任务比较分析函数
compare_tasks <- function(gam_results, successful_metrics) {
  
  cat("\n=== Inter-task Comparison Analysis ===\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # Group by task
  stroop_metrics <- successful_metrics[grepl("Stroop", successful_metrics)]
  gonogo_metrics <- successful_metrics[grepl("GoNoGo", successful_metrics)]
  
  # Calculate task averages
  calculate_task_average <- function(metrics, results_list, var_name) {
    if (length(metrics) == 0) return(NA)
    values <- sapply(metrics, function(m) results_list[[m]][[var_name]])
    return(mean(values, na.rm = TRUE))
  }
  
  cat("Stroop task (n =", length(stroop_metrics), "):\n")
  cat("  Average R²: ", round(calculate_task_average(stroop_metrics, gam_results, "deviance_explained") * 100, 1), "%\n")
  cat("  Average EDF: ", round(calculate_task_average(stroop_metrics, gam_results, "edf"), 2), "\n")
  cat("  Average AIC: ", round(calculate_task_average(stroop_metrics, gam_results, "aic"), 2), "\n")
  cat("  Average F-value: ", round(calculate_task_average(stroop_metrics, gam_results, "smooth_f_value"), 2), "\n")
  cat("  Average p-value: ", round(calculate_task_average(stroop_metrics, gam_results, "smooth_p_value"), 4), "\n")
  
  cat("\nGoNoGo task (n =", length(gonogo_metrics), "):\n")
  cat("  Average R²: ", round(calculate_task_average(gonogo_metrics, gam_results, "deviance_explained") * 100, 1), "%\n")
  cat("  Average EDF: ", round(calculate_task_average(gonogo_metrics, gam_results, "edf"), 2), "\n")
  cat("  Average AIC: ", round(calculate_task_average(gonogo_metrics, gam_results, "aic"), 2), "\n")
  cat("  Average F-value: ", round(calculate_task_average(gonogo_metrics, gam_results, "smooth_f_value"), 2), "\n")
  cat("  Average p-value: ", round(calculate_task_average(gonogo_metrics, gam_results, "smooth_p_value"), 4), "\n")
  
  # Core inhibitory control metrics comparison
  core_inhibition <- successful_metrics[grepl("Stroop_Incongruent|GoNoGo_Nogo_ACC|GoNoGo_Gonogo_d'|Gonogo_d", successful_metrics)]
  
  if (length(core_inhibition) > 0) {
    cat("\nCore inhibitory control metrics (n =", length(core_inhibition), "):\n")
    cat("  Average R²: ", round(calculate_task_average(core_inhibition, gam_results, "deviance_explained") * 100, 1), "%\n")
    cat("  Average EDF: ", round(calculate_task_average(core_inhibition, gam_results, "edf"), 2), "\n")
    cat("  Average AIC: ", round(calculate_task_average(core_inhibition, gam_results, "aic"), 2), "\n")
    cat("  Average F-value: ", round(calculate_task_average(core_inhibition, gam_results, "smooth_f_value"), 2), "\n")
    cat("  Average p-value: ", round(calculate_task_average(core_inhibition, gam_results, "smooth_p_value"), 4), "\n")
  }
  
  # Significance statistics
  cat("\n=== Smooth Term Significance Statistics ===\n")
  significant_001 <- sum(sapply(successful_metrics, function(m) gam_results[[m]]$smooth_p_value < 0.001))
  significant_01 <- sum(sapply(successful_metrics, function(m) gam_results[[m]]$smooth_p_value < 0.01))
  significant_05 <- sum(sapply(successful_metrics, function(m) gam_results[[m]]$smooth_p_value < 0.05))
  significant_1 <- sum(sapply(successful_metrics, function(m) gam_results[[m]]$smooth_p_value < 0.1))
  
  cat("Metrics with p < 0.001: ", significant_001, " (", round(significant_001/length(successful_metrics)*100, 1), "%)\n")
  cat("Metrics with p < 0.01: ", significant_01, " (", round(significant_01/length(successful_metrics)*100, 1), "%)\n")
  cat("Metrics with p < 0.05: ", significant_05, " (", round(significant_05/length(successful_metrics)*100, 1), "%)\n")
  cat("Metrics with p < 0.1: ", significant_1, " (", round(significant_1/length(successful_metrics)*100, 1), "%)\n")
}

# 保存核心指标图形函数
save_core_plots <- function(core_plots, output_dir = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/output") {
  if (is.null(core_plots) || length(core_plots) == 0) {
    cat("No core metric plots to save\n")
    return()
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  core_dir <- file.path(output_dir, "core_metrics_plots")
  if (!dir.exists(core_dir)) {
    dir.create(core_dir, recursive = TRUE)
  }
  
  for (metric_name in names(core_plots)) {
    tryCatch({
      filename <- file.path(core_dir, paste0(metric_name, "_individual_fit.png"))
      ggsave(filename, core_plots[[metric_name]],
             width = 8, height = 5, dpi = 800)
      cat("Saved core metric plot:", filename, "\n")
    }, error = function(e) {
      cat("Failed to save", metric_name, "core metric plot:", e$message, "\n")
    })
  }
}

# 创建平滑项显著性汇总表函数
create_significance_summary <- function(gam_results, successful_metrics) {
  
  cat("\n=== Creating Smooth Term Significance Summary ===\n")
  
  sig_summary <- data.frame()
  
  for (metric in successful_metrics) {
    results <- gam_results[[metric]]
    
    # Determine significance level
    p_val <- results$smooth_p_value
    if (p_val < 0.001) {
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
      sig_text <- "Not significant"
    }
    
    row_data <- data.frame(
      Metric = metric,
      Task = ifelse(grepl("Stroop", metric), "Stroop", "GoNoGo"),
      Variable = gsub("^(Stroop_|GoNoGo_)", "", metric),
      F_Value = round(results$smooth_f_value, 3),
      P_Value = round(results$smooth_p_value, 6),
      Ref_DF = round(results$smooth_ref_df, 2),
      EDF = round(results$edf, 2),
      Significance_Level = sig_level,
      Significance_Text = sig_text,
      R_squared = round(results$deviance_explained * 100, 1),
      AIC = round(results$aic, 2),
      stringsAsFactors = FALSE
    )
    
    sig_summary <- rbind(sig_summary, row_data)
  }
  
  # Sort by task and significance
  sig_summary <- sig_summary[order(sig_summary$Task, sig_summary$P_Value), ]
  
  return(sig_summary)
}

# 创建诊断图函数
create_diagnostic_plots <- function(gam_results, successful_metrics) {
  
  # Placeholder for diagnostic plots
  # Can be expanded to include residual plots, QQ plots, etc.
  diagnostic_plots <- list()
  
  # Simple return empty list, can be expanded as needed
  return(diagnostic_plots)
}

# 保存诊断图函数
save_diagnostic_plots <- function(diagnostic_plots, output_dir) {
  
  # Placeholder for saving diagnostic plots
  # Simple return, can be expanded as needed
  return()
}

# 执行完整分析并生成报告
run_complete_analysis <- function(stroop_path = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Stroop-baseline.xlsx",
                                  gonogo_path = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Gonogo-baseline.xlsx",
                                  output_dir = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/output") {
  
  cat("=== Starting Complete GAM Analysis Pipeline ===\n")
  
  # 1. Run main analysis
  results <- main_gam_analysis(stroop_path, gonogo_path)
  
  # 2. Run model diagnostics
  run_model_diagnostics(results$gam_results, results$successful_metrics)
  
  # 3. Task comparison analysis
  compare_tasks(results$gam_results, results$successful_metrics)
  
  # 4. Create smooth term significance summary
  sig_summary <- create_significance_summary(results$gam_results, results$successful_metrics)
  
  # 5. Create diagnostic plots
  diagnostic_plots <- create_diagnostic_plots(results$gam_results, results$successful_metrics)
  
  # 6. Save diagnostic plots
  save_diagnostic_plots(diagnostic_plots, output_dir)
  
  # 6.5. Save core metric plots
  if (!is.null(results$plots$core_plots)) {
    save_core_plots(results$plots$core_plots, output_dir)
  }
  
  # 7. Save significance summary
  tryCatch({
    write.csv(sig_summary, file.path(output_dir, "smooth_significance_summary.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat("\nSignificance summary saved to:", file.path(output_dir, "smooth_significance_summary.csv"), "\n")
  }, error = function(e) {
    cat("Failed to save significance summary:", e$message, "\n")
  })
  
  # 8. Print significance summary
  cat("\n=== Smooth Term Significance Summary ===\n")
  cat(rep("=", 80), "\n", sep = "")
  print(sig_summary[, c("Metric", "F_Value", "P_Value", "Significance_Level",
                        "EDF", "R_squared")], row.names = FALSE)
  
  cat("\n=== Complete Analysis Pipeline Finished ===\n")
  
  return(list(
    main_results = results,
    significance_summary = sig_summary,
    diagnostic_plots = diagnostic_plots
  ))
}

# 执行完整分析
cat("R GAM Analysis Script Ready\n")
cat("Usage:\n")
cat("1. Run main analysis: results <- main_gam_analysis()\n")
cat("2. Run model diagnostics: run_model_diagnostics(results$gam_results, results$successful_metrics)\n")
cat("3. Task comparison: compare_tasks(results$gam_results, results$successful_metrics)\n")
cat("4. Complete analysis: complete_results <- run_complete_analysis()\n")

# 检查文件是否存在
stroop_file <- "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Stroop-baseline.xlsx"
gonogo_file <- "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Gonogo-baseline.xlsx"

cat("\n=== 文件检查 ===\n")
cat("Stroop文件路径:", stroop_file, "\n")
cat("Stroop文件存在:", file.exists(stroop_file), "\n")
cat("GoNoGo文件路径:", gonogo_file, "\n")
cat("GoNoGo文件存在:", file.exists(gonogo_file), "\n")

if (file.exists(stroop_file) && file.exists(gonogo_file)) {
  cat("\n所有文件都存在，开始分析...\n")
  # Run analysis (uncomment to execute)
  complete_results <- run_complete_analysis()
} else {
  cat("\n错误：某些文件不存在，请检查文件路径\n")
  cat("当前工作目录:", getwd(), "\n")
  cat("请使用以下命令之一:\n")
  cat("1. 设置正确的工作目录: setwd('您的文件夹路径')\n")
  cat("2. 或者使用完整路径运行: complete_results <- run_complete_analysis('完整的Stroop路径', '完整的GoNoGo路径')\n")
}

# ========== 附加功能函数 ==========

# 创建热力图可视化函数
create_heatmap_visualization <- function(gam_results, successful_metrics, output_dir = NULL) {
  
  cat("\n=== 创建发展轨迹热力图 ===\n")
  
  if (length(successful_metrics) == 0) {
    cat("没有成功的指标用于创建热力图\n")
    return(NULL)
  }
  
  # 准备热力图数据
  age_seq <- seq(6, 18, by = 0.1)
  heatmap_data <- data.frame(age = age_seq)
  
  # 收集所有指标的预测值和导数
  for (metric in successful_metrics) {
    results <- gam_results[[metric]]
    if (!is.null(results)) {
      # 插值到统一网格
      pred_interp <- approx(results$age, results$predictions, xout = age_seq, rule = 2)$y
      deriv_interp <- approx(results$age, results$derivatives, xout = age_seq, rule = 2)$y
      sig_interp <- approx(results$age, as.numeric(results$significant), xout = age_seq, rule = 2)$y
      
      heatmap_data[[paste0(metric, "_pred")]] <- pred_interp
      heatmap_data[[paste0(metric, "_deriv")]] <- deriv_interp
      heatmap_data[[paste0(metric, "_sig")]] <- sig_interp >= 0.5
    }
  }
  
  # 转换为长格式用于ggplot
  pred_cols <- grep("_pred$", names(heatmap_data), value = TRUE)
  deriv_cols <- grep("_deriv$", names(heatmap_data), value = TRUE)
  sig_cols <- grep("_sig$", names(heatmap_data), value = TRUE)
  
  # 预测值热力图数据
  pred_long <- heatmap_data %>%
    select(age, all_of(pred_cols)) %>%
    gather(key = "metric", value = "prediction", -age) %>%
    mutate(metric = gsub("_pred$", "", metric))
  
  # 导数热力图数据
  deriv_long <- heatmap_data %>%
    select(age, all_of(deriv_cols)) %>%
    gather(key = "metric", value = "derivative", -age) %>%
    mutate(metric = gsub("_deriv$", "", metric))
  
  # 显著性数据
  sig_long <- heatmap_data %>%
    select(age, all_of(sig_cols)) %>%
    gather(key = "metric", value = "significant", -age) %>%
    mutate(metric = gsub("_sig$", "", metric))
  
  # 创建预测值热力图
  p_pred_heatmap <- ggplot(pred_long, aes(x = age, y = metric, fill = prediction)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, name = "Z-score") +
    labs(x = "Age (years)", y = "Metrics", 
         title = "Development Trajectory Heatmap: Predicted Values") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 创建导数热力图
  p_deriv_heatmap <- ggplot(deriv_long, aes(x = age, y = metric, fill = derivative)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = 0, name = "Change Rate") +
    labs(x = "Age (years)", y = "Metrics", 
         title = "Development Rate Heatmap: Derivatives") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 创建显著性热力图
  p_sig_heatmap <- ggplot(sig_long, aes(x = age, y = metric, fill = significant)) +
    geom_tile() +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "darkred"), 
                      name = "Significant") +
    labs(x = "Age (years)", y = "Metrics", 
         title = "Significant Development Periods") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 保存热力图
  if (!is.null(output_dir)) {
    tryCatch({
      ggsave(file.path(output_dir, "prediction_heatmap.png"), 
             p_pred_heatmap, width = 12, height = 8, dpi = 300)
      ggsave(file.path(output_dir, "derivative_heatmap.png"), 
             p_deriv_heatmap, width = 12, height = 8, dpi = 300)
      ggsave(file.path(output_dir, "significance_heatmap.png"), 
             p_sig_heatmap, width = 12, height = 8, dpi = 300)
      cat("热力图已保存到:", output_dir, "\n")
    }, error = function(e) {
      cat("保存热力图失败:", e$message, "\n")
    })
  }
  
  return(list(
    prediction_heatmap = p_pred_heatmap,
    derivative_heatmap = p_deriv_heatmap,
    significance_heatmap = p_sig_heatmap,
    heatmap_data = heatmap_data
  ))
}

# 创建发展阶段总结报告函数
create_development_summary_report <- function(gam_results, successful_metrics, 
                                              inverted_u_stats, stage_stats, 
                                              completion_stats, output_dir = NULL) {
  
  cat("\n=== 生成发展阶段总结报告 ===\n")
  
  # 准备报告数据
  report_data <- data.frame()
  
  for (metric in successful_metrics) {
    # 基本信息
    results <- gam_results[[metric]]
    u_stats <- inverted_u_stats[[metric]]
    stage_info <- stage_stats[[metric]]
    completion_info <- completion_stats[[metric]]
    
    # 任务和变量信息
    task_type <- ifelse(grepl("Stroop", metric), "Stroop", "GoNoGo")
    variable_name <- gsub("^(Stroop_|GoNoGo_)", "", metric)
    
    # 显著发展期信息
    sig_periods <- results$significant
    if (any(sig_periods)) {
      sig_ages <- results$age[sig_periods]
      sig_start <- min(sig_ages)
      sig_end <- max(sig_ages)
      sig_duration <- sig_end - sig_start
      sig_coverage <- sum(sig_periods) / length(sig_periods) * 100
    } else {
      sig_start <- NA
      sig_end <- NA
      sig_duration <- 0
      sig_coverage <- 0
    }
    
    # 峰值变化信息
    max_change_idx <- which.max(abs(results$derivatives))
    peak_change_age <- results$age[max_change_idx]
    peak_change_rate <- results$derivatives[max_change_idx]
    
    # 创建行数据
    row_data <- data.frame(
      Metric = metric,
      Task = task_type,
      Variable = variable_name,
      
      # 模型拟合质量
      R_Squared = round(results$deviance_explained * 100, 1),
      AIC = round(results$aic, 1),
      F_Value = round(results$smooth_f_value, 2),
      P_Value = round(results$smooth_p_value, 4),
      
      # 显著发展期
      Sig_Start_Age = ifelse(is.na(sig_start), "无", round(sig_start, 1)),
      Sig_End_Age = ifelse(is.na(sig_end), "无", round(sig_end, 1)),
      Sig_Duration_Years = round(sig_duration, 1),
      Sig_Coverage_Percent = round(sig_coverage, 1),
      
      # 峰值变化
      Peak_Change_Age = round(peak_change_age, 1),
      Peak_Change_Rate = round(peak_change_rate, 3),
      
      # 倒U型模式
      Has_Inverted_U = ifelse(u_stats$has_inverted_u, "是", "否"),
      U_Peak_Age = ifelse(is.na(u_stats$peak_age), "无", round(u_stats$peak_age, 1)),
      
      # 发展阶段持续时间
      Early_Stage_Duration = round(stage_info$early_rapid_6_8_duration, 1),
      Middle_Stage_Duration = round(stage_info$middle_sustained_9_12_duration, 1),
      Late_Stage_Duration = round(stage_info$late_plateau_13_15_duration, 1),
      Final_Stage_Duration = round(stage_info$final_stable_16_18_duration, 1),
      
      # 12岁前完成度
      Completion_12 = round(completion_info$completion_12, 1),
      Completion_15 = round(completion_info$completion_15, 1),
      
      stringsAsFactors = FALSE
    )
    
    report_data <- rbind(report_data, row_data)
  }
  
  # 保存报告
  if (!is.null(output_dir)) {
    tryCatch({
      write.csv(report_data, file.path(output_dir, "development_summary_report.csv"),
                row.names = FALSE, fileEncoding = "UTF-8")
      cat("发展总结报告已保存到:", file.path(output_dir, "development_summary_report.csv"), "\n")
    }, error = function(e) {
      cat("保存发展总结报告失败:", e$message, "\n")
    })
  }
  
  # 打印关键发现
  cat("\n=== 关键发现总结 ===\n")
  cat(rep("=", 50), "\n", sep = "")
  
  # 按任务类型分组统计
  stroop_metrics <- report_data[report_data$Task == "Stroop", ]
  gonogo_metrics <- report_data[report_data$Task == "GoNoGo", ]
  
  if (nrow(stroop_metrics) > 0) {
    cat("\nStroop任务 (", nrow(stroop_metrics), "个指标):\n")
    cat("  平均R²: ", round(mean(as.numeric(stroop_metrics$R_Squared), na.rm = TRUE), 1), "%\n")
    cat("  有显著发展期的指标: ", sum(stroop_metrics$Sig_Duration_Years > 0), "个\n")
    cat("  倒U型发展模式: ", sum(stroop_metrics$Has_Inverted_U == "是"), "个\n")
    cat("  平均显著发展持续时间: ", round(mean(as.numeric(stroop_metrics$Sig_Duration_Years), na.rm = TRUE), 1), "年\n")
  }
  
  if (nrow(gonogo_metrics) > 0) {
    cat("\nGoNoGo任务 (", nrow(gonogo_metrics), "个指标):\n")
    cat("  平均R²: ", round(mean(as.numeric(gonogo_metrics$R_Squared), na.rm = TRUE), 1), "%\n")
    cat("  有显著发展期的指标: ", sum(gonogo_metrics$Sig_Duration_Years > 0), "个\n")
    cat("  倒U型发展模式: ", sum(gonogo_metrics$Has_Inverted_U == "是"), "个\n")
    cat("  平均显著发展持续时间: ", round(mean(as.numeric(gonogo_metrics$Sig_Duration_Years), na.rm = TRUE), 1), "年\n")
  }
  
  # 核心抑制控制指标分析
  core_inhibition <- report_data[grepl("Incongruent|Nogo_ACC|Gonogo_d|d'", report_data$Variable), ]
  if (nrow(core_inhibition) > 0) {
    cat("\n核心抑制控制指标 (", nrow(core_inhibition), "个指标):\n")
    cat("  平均R²: ", round(mean(as.numeric(core_inhibition$R_Squared), na.rm = TRUE), 1), "%\n")
    cat("  有显著发展期的指标: ", sum(core_inhibition$Sig_Duration_Years > 0), "个\n")
    cat("  12岁前平均完成度: ", round(mean(as.numeric(core_inhibition$Completion_12), na.rm = TRUE), 1), "%\n")
    cat("  15岁前平均完成度: ", round(mean(as.numeric(core_inhibition$Completion_15), na.rm = TRUE), 1), "%\n")
  }
  
  return(report_data)
}

# 扩展的完整分析函数
run_extended_complete_analysis <- function(stroop_path = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Stroop-baseline.xlsx",
                                           gonogo_path = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/Gonogo-baseline.xlsx",
                                           output_dir = "/Users/zhangtongyi/Desktop/博士论文/1.博士论文计算建模研究/2.一般分析/2.Gonogo和Stroop比-跨站点/output") {
  
  cat("=== 开始扩展完整GAM分析流程 ===\n")
  
  # 1. 运行基础分析
  basic_results <- run_complete_analysis(stroop_path, gonogo_path, output_dir)
  
  # 2. 创建热力图可视化
  heatmap_results <- create_heatmap_visualization(
    basic_results$main_results$gam_results, 
    basic_results$main_results$successful_metrics, 
    output_dir
  )
  
  # 3. 生成发展阶段总结报告
  summary_report <- create_development_summary_report(
    basic_results$main_results$gam_results,
    basic_results$main_results$successful_metrics,
    basic_results$main_results$inverted_u_stats,
    basic_results$main_results$stage_stats,
    basic_results$main_results$completion_stats,
    output_dir
  )
  
  cat("\n=== 扩展完整分析流程完成 ===\n")
  cat("所有结果和图形已保存到:", output_dir, "\n")
  
  return(list(
    basic_results = basic_results,
    heatmap_results = heatmap_results,
    summary_report = summary_report
  ))
}

# ========== 快速启动命令 ==========

cat("\n" , rep("=", 60), "\n", sep = "")
cat("GAM 分析脚本已完全加载！\n")
cat(rep("=", 60), "\n", sep = "")
cat("\n快速启动命令:\n")
cat("1. 基础分析: \n")
cat("   results <- main_gam_analysis()\n\n")
cat("2. 完整分析: \n")
cat("   complete_results <- run_complete_analysis()\n\n")
cat("3. 扩展完整分析（包括热力图）: \n")
cat("   extended_results <- run_extended_complete_analysis()\n\n")
cat("4. 单独创建热力图: \n")
cat("   heatmaps <- create_heatmap_visualization(results$gam_results, results$successful_metrics)\n\n")
cat("5. 模型诊断: \n")
cat("   run_model_diagnostics(results$gam_results, results$successful_metrics)\n\n")
cat(rep("=", 60), "\n", sep = "")
