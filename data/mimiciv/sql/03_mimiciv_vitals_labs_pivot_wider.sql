CREATE OR REPLACE TABLE `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_24h_features_wide_3` AS
SELECT
  stay_id,
  mortality,
  gender,
  anchor_age,

  -- =========================
  -- HEART RATE
  -- =========================
  MAX(IF(signal = 'heart_rate', value_min, NULL)) AS heart_rate_min,
  MAX(IF(signal = 'heart_rate', value_max, NULL)) AS heart_rate_max,
  MAX(IF(signal = 'heart_rate', value_mean, NULL)) AS heart_rate_mean,
  MAX(IF(signal = 'heart_rate', value_median, NULL)) AS heart_rate_median,
  MAX(IF(signal = 'heart_rate', value_iqr, NULL)) AS heart_rate_iqr,
  MAX(IF(signal = 'heart_rate', measurement_count, NULL)) AS heart_rate_count,
  MAX(IF(signal = 'heart_rate', missing_indicator, NULL)) AS heart_rate_missing,
  MAX(IF(signal = 'heart_rate', early_mean, NULL)) AS heart_rate_early_mean,
  MAX(IF(signal = 'heart_rate', late_mean, NULL)) AS heart_rate_late_mean,
  MAX(IF(signal = 'heart_rate', count_below_safe_min, NULL)) AS heart_rate_low_count,
  MAX(IF(signal = 'heart_rate', count_above_safe_max, NULL)) AS heart_rate_high_count,
  MAX(IF(signal = 'heart_rate', bin_low_count, NULL)) AS heart_rate_bin_low,
  MAX(IF(signal = 'heart_rate', bin_normal_count, NULL)) AS heart_rate_bin_normal,
  MAX(IF(signal = 'heart_rate', bin_high_count, NULL)) AS heart_rate_bin_high,

  -- =========================
  -- MBP
  -- =========================
  MAX(IF(signal = 'mbp', value_min, NULL)) AS mbp_min,
  MAX(IF(signal = 'mbp', value_max, NULL)) AS mbp_max,
  MAX(IF(signal = 'mbp', value_mean, NULL)) AS mbp_mean,
  MAX(IF(signal = 'mbp', value_median, NULL)) AS mbp_median,
  MAX(IF(signal = 'mbp', value_iqr, NULL)) AS mbp_iqr,
  MAX(IF(signal = 'mbp', measurement_count, NULL)) AS mbp_count,
  MAX(IF(signal = 'mbp', missing_indicator, NULL)) AS mbp_missing,
  MAX(IF(signal = 'mbp', early_mean, NULL)) AS mbp_early_mean,
  MAX(IF(signal = 'mbp', late_mean, NULL)) AS mbp_late_mean,
  MAX(IF(signal = 'mbp', count_below_safe_min, NULL)) AS mbp_low_count,
  MAX(IF(signal = 'mbp', count_above_safe_max, NULL)) AS mbp_high_count,
  MAX(IF(signal = 'mbp', bin_low_count, NULL)) AS mbp_bin_low,
  MAX(IF(signal = 'mbp', bin_normal_count, NULL)) AS mbp_bin_normal,
  MAX(IF(signal = 'mbp', bin_high_count, NULL)) AS mbp_bin_high,

  -- =========================
  -- SBP
  -- =========================
  MAX(IF(signal = 'sbp', value_min, NULL)) AS sbp_min,
  MAX(IF(signal = 'sbp', value_max, NULL)) AS sbp_max,
  MAX(IF(signal = 'sbp', value_mean, NULL)) AS sbp_mean,
  MAX(IF(signal = 'sbp', value_median, NULL)) AS sbp_median,
  MAX(IF(signal = 'sbp', value_iqr, NULL)) AS sbp_iqr,
  MAX(IF(signal = 'sbp', measurement_count, NULL)) AS sbp_count,
  MAX(IF(signal = 'sbp', missing_indicator, NULL)) AS sbp_missing,
  MAX(IF(signal = 'sbp', early_mean, NULL)) AS sbp_early_mean,
  MAX(IF(signal = 'sbp', late_mean, NULL)) AS sbp_late_mean,
  MAX(IF(signal = 'sbp', count_below_safe_min, NULL)) AS sbp_low_count,
  MAX(IF(signal = 'sbp', count_above_safe_max, NULL)) AS sbp_high_count,
  MAX(IF(signal = 'sbp', bin_low_count, NULL)) AS sbp_bin_low,
  MAX(IF(signal = 'sbp', bin_normal_count, NULL)) AS sbp_bin_normal,
  MAX(IF(signal = 'sbp', bin_high_count, NULL)) AS sbp_bin_high,

  -- =========================
  -- DBP
  -- =========================
  MAX(IF(signal = 'dbp', value_min, NULL)) AS dbp_min,
  MAX(IF(signal = 'dbp', value_max, NULL)) AS dbp_max,
  MAX(IF(signal = 'dbp', value_mean, NULL)) AS dbp_mean,
  MAX(IF(signal = 'dbp', value_median, NULL)) AS dbp_median,
  MAX(IF(signal = 'dbp', value_iqr, NULL)) AS dbp_iqr,
  MAX(IF(signal = 'dbp', measurement_count, NULL)) AS dbp_count,
  MAX(IF(signal = 'dbp', missing_indicator, NULL)) AS dbp_missing,
  MAX(IF(signal = 'dbp', early_mean, NULL)) AS dbp_early_mean,
  MAX(IF(signal = 'dbp', late_mean, NULL)) AS dbp_late_mean,
  MAX(IF(signal = 'dbp', count_below_safe_min, NULL)) AS dbp_low_count,
  MAX(IF(signal = 'dbp', count_above_safe_max, NULL)) AS dbp_high_count,
  MAX(IF(signal = 'dbp', bin_low_count, NULL)) AS dbp_bin_low,
  MAX(IF(signal = 'dbp', bin_normal_count, NULL)) AS dbp_bin_normal,
  MAX(IF(signal = 'dbp', bin_high_count, NULL)) AS dbp_bin_high,

  -- =========================
  -- RESPIRATORY RATE
  -- =========================
  MAX(IF(signal = 'resp_rate', value_min, NULL)) AS resp_rate_min,
  MAX(IF(signal = 'resp_rate', value_max, NULL)) AS resp_rate_max,
  MAX(IF(signal = 'resp_rate', value_mean, NULL)) AS resp_rate_mean,
  MAX(IF(signal = 'resp_rate', value_median, NULL)) AS resp_rate_median,
  MAX(IF(signal = 'resp_rate', value_iqr, NULL)) AS resp_rate_iqr,
  MAX(IF(signal = 'resp_rate', measurement_count, NULL)) AS resp_rate_count,
  MAX(IF(signal = 'resp_rate', missing_indicator, NULL)) AS resp_rate_missing,
  MAX(IF(signal = 'resp_rate', early_mean, NULL)) AS resp_rate_early_mean,
  MAX(IF(signal = 'resp_rate', late_mean, NULL)) AS resp_rate_late_mean,
  MAX(IF(signal = 'resp_rate', count_below_safe_min, NULL)) AS resp_rate_low_count,
  MAX(IF(signal = 'resp_rate', count_above_safe_max, NULL)) AS resp_rate_high_count,
  MAX(IF(signal = 'resp_rate', bin_low_count, NULL)) AS resp_rate_bin_low,
  MAX(IF(signal = 'resp_rate', bin_normal_count, NULL)) AS resp_rate_bin_normal,
  MAX(IF(signal = 'resp_rate', bin_high_count, NULL)) AS resp_rate_bin_high,

  -- =========================
  -- SPO2
  -- =========================
  MAX(IF(signal = 'spo2', value_min, NULL)) AS spo2_min,
  MAX(IF(signal = 'spo2', value_max, NULL)) AS spo2_max,
  MAX(IF(signal = 'spo2', value_mean, NULL)) AS spo2_mean,
  MAX(IF(signal = 'spo2', value_median, NULL)) AS spo2_median,
  MAX(IF(signal = 'spo2', value_iqr, NULL)) AS spo2_iqr,
  MAX(IF(signal = 'spo2', measurement_count, NULL)) AS spo2_count,
  MAX(IF(signal = 'spo2', missing_indicator, NULL)) AS spo2_missing,
  MAX(IF(signal = 'spo2', early_mean, NULL)) AS spo2_early_mean,
  MAX(IF(signal = 'spo2', late_mean, NULL)) AS spo2_late_mean,
  MAX(IF(signal = 'spo2', count_below_safe_min, NULL)) AS spo2_low_count,
  MAX(IF(signal = 'spo2', count_above_safe_max, NULL)) AS spo2_high_count,
  MAX(IF(signal = 'spo2', bin_low_count, NULL)) AS spo2_bin_low,
  MAX(IF(signal = 'spo2', bin_normal_count, NULL)) AS spo2_bin_normal,
  MAX(IF(signal = 'spo2', bin_high_count, NULL)) AS spo2_bin_high,

  -- =========================
  -- TEMPERATURE
  -- =========================
  MAX(IF(signal = 'temperature', value_min, NULL)) AS temperature_min,
  MAX(IF(signal = 'temperature', value_max, NULL)) AS temperature_max,
  MAX(IF(signal = 'temperature', value_mean, NULL)) AS temperature_mean,
  MAX(IF(signal = 'temperature', value_median, NULL)) AS temperature_median,
  MAX(IF(signal = 'temperature', value_iqr, NULL)) AS temperature_iqr,
  MAX(IF(signal = 'temperature', measurement_count, NULL)) AS temperature_count,
  MAX(IF(signal = 'temperature', missing_indicator, NULL)) AS temperature_missing,
  MAX(IF(signal = 'temperature', early_mean, NULL)) AS temperature_early_mean,
  MAX(IF(signal = 'temperature', late_mean, NULL)) AS temperature_late_mean,
  MAX(IF(signal = 'temperature', count_below_safe_min, NULL)) AS temperature_low_count,
  MAX(IF(signal = 'temperature', count_above_safe_max, NULL)) AS temperature_high_count,
  MAX(IF(signal = 'temperature', bin_low_count, NULL)) AS temperature_bin_low,
  MAX(IF(signal = 'temperature', bin_normal_count, NULL)) AS temperature_bin_normal,
  MAX(IF(signal = 'temperature', bin_high_count, NULL)) AS temperature_bin_high,

  -- =========================
  -- CREATININE
  -- =========================
  MAX(IF(signal = 'creatinine', value_min, NULL)) AS creatinine_min,
  MAX(IF(signal = 'creatinine', value_max, NULL)) AS creatinine_max,
  MAX(IF(signal = 'creatinine', value_mean, NULL)) AS creatinine_mean,
  MAX(IF(signal = 'creatinine', value_median, NULL)) AS creatinine_median,
  MAX(IF(signal = 'creatinine', value_iqr, NULL)) AS creatinine_iqr,
  MAX(IF(signal = 'creatinine', measurement_count, NULL)) AS creatinine_count,
  MAX(IF(signal = 'creatinine', missing_indicator, NULL)) AS creatinine_missing,
  MAX(IF(signal = 'creatinine', early_mean, NULL)) AS creatinine_early_mean,
  MAX(IF(signal = 'creatinine', late_mean, NULL)) AS creatinine_late_mean,
  MAX(IF(signal = 'creatinine', count_below_safe_min, NULL)) AS creatinine_low_count,
  MAX(IF(signal = 'creatinine', count_above_safe_max, NULL)) AS creatinine_high_count,
  MAX(IF(signal = 'creatinine', bin_low_count, NULL)) AS creatinine_bin_low,
  MAX(IF(signal = 'creatinine', bin_normal_count, NULL)) AS creatinine_bin_normal,
  MAX(IF(signal = 'creatinine', bin_high_count, NULL)) AS creatinine_bin_high,

  -- =========================
  -- SODIUM
  -- =========================
  MAX(IF(signal = 'sodium', value_min, NULL)) AS sodium_min,
  MAX(IF(signal = 'sodium', value_max, NULL)) AS sodium_max,
  MAX(IF(signal = 'sodium', value_mean, NULL)) AS sodium_mean,
  MAX(IF(signal = 'sodium', value_median, NULL)) AS sodium_median,
  MAX(IF(signal = 'sodium', value_iqr, NULL)) AS sodium_iqr,
  MAX(IF(signal = 'sodium', measurement_count, NULL)) AS sodium_count,
  MAX(IF(signal = 'sodium', missing_indicator, NULL)) AS sodium_missing,
  MAX(IF(signal = 'sodium', early_mean, NULL)) AS sodium_early_mean,
  MAX(IF(signal = 'sodium', late_mean, NULL)) AS sodium_late_mean,
  MAX(IF(signal = 'sodium', count_below_safe_min, NULL)) AS sodium_low_count,
  MAX(IF(signal = 'sodium', count_above_safe_max, NULL)) AS sodium_high_count,
  MAX(IF(signal = 'sodium', bin_low_count, NULL)) AS sodium_bin_low,
  MAX(IF(signal = 'sodium', bin_normal_count, NULL)) AS sodium_bin_normal,
  MAX(IF(signal = 'sodium', bin_high_count, NULL)) AS sodium_bin_high,

  -- =========================
  -- POTASSIUM
  -- =========================
  MAX(IF(signal = 'potassium', value_min, NULL)) AS potassium_min,
  MAX(IF(signal = 'potassium', value_max, NULL)) AS potassium_max,
  MAX(IF(signal = 'potassium', value_mean, NULL)) AS potassium_mean,
  MAX(IF(signal = 'potassium', value_median, NULL)) AS potassium_median,
  MAX(IF(signal = 'potassium', value_iqr, NULL)) AS potassium_iqr,
  MAX(IF(signal = 'potassium', measurement_count, NULL)) AS potassium_count,
  MAX(IF(signal = 'potassium', missing_indicator, NULL)) AS potassium_missing,
  MAX(IF(signal = 'potassium', early_mean, NULL)) AS potassium_early_mean,
  MAX(IF(signal = 'potassium', late_mean, NULL)) AS potassium_late_mean,
  MAX(IF(signal = 'potassium', count_below_safe_min, NULL)) AS potassium_low_count,
  MAX(IF(signal = 'potassium', count_above_safe_max, NULL)) AS potassium_high_count,
  MAX(IF(signal = 'potassium', bin_low_count, NULL)) AS potassium_bin_low,
  MAX(IF(signal = 'potassium', bin_normal_count, NULL)) AS potassium_bin_normal,
  MAX(IF(signal = 'potassium', bin_high_count, NULL)) AS potassium_bin_high,

  -- =========================
  -- BICARBONATE
  -- =========================
  MAX(IF(signal = 'bicarbonate', value_min, NULL)) AS bicarbonate_min,
  MAX(IF(signal = 'bicarbonate', value_max, NULL)) AS bicarbonate_max,
  MAX(IF(signal = 'bicarbonate', value_mean, NULL)) AS bicarbonate_mean,
  MAX(IF(signal = 'bicarbonate', value_median, NULL)) AS bicarbonate_median,
  MAX(IF(signal = 'bicarbonate', value_iqr, NULL)) AS bicarbonate_iqr,
  MAX(IF(signal = 'bicarbonate', measurement_count, NULL)) AS bicarbonate_count,
  MAX(IF(signal = 'bicarbonate', missing_indicator, NULL)) AS bicarbonate_missing,
  MAX(IF(signal = 'bicarbonate', early_mean, NULL)) AS bicarbonate_early_mean,
  MAX(IF(signal = 'bicarbonate', late_mean, NULL)) AS bicarbonate_late_mean,
  MAX(IF(signal = 'bicarbonate', count_below_safe_min, NULL)) AS bicarbonate_low_count,
  MAX(IF(signal = 'bicarbonate', count_above_safe_max, NULL)) AS bicarbonate_high_count,
  MAX(IF(signal = 'bicarbonate', bin_low_count, NULL)) AS bicarbonate_bin_low,
  MAX(IF(signal = 'bicarbonate', bin_normal_count, NULL)) AS bicarbonate_bin_normal,
  MAX(IF(signal = 'bicarbonate', bin_high_count, NULL)) AS bicarbonate_bin_high,

  -- =========================
  -- HEMOGLOBIN
  -- =========================
  MAX(IF(signal = 'hemoglobin', value_min, NULL)) AS hemoglobin_min,
  MAX(IF(signal = 'hemoglobin', value_max, NULL)) AS hemoglobin_max,
  MAX(IF(signal = 'hemoglobin', value_mean, NULL)) AS hemoglobin_mean,
  MAX(IF(signal = 'hemoglobin', value_median, NULL)) AS hemoglobin_median,
  MAX(IF(signal = 'hemoglobin', value_iqr, NULL)) AS hemoglobin_iqr,
  MAX(IF(signal = 'hemoglobin', measurement_count, NULL)) AS hemoglobin_count,
  MAX(IF(signal = 'hemoglobin', missing_indicator, NULL)) AS hemoglobin_missing,
  MAX(IF(signal = 'hemoglobin', early_mean, NULL)) AS hemoglobin_early_mean,
  MAX(IF(signal = 'hemoglobin', late_mean, NULL)) AS hemoglobin_late_mean,
  MAX(IF(signal = 'hemoglobin', count_below_safe_min, NULL)) AS hemoglobin_low_count,
  MAX(IF(signal = 'hemoglobin', count_above_safe_max, NULL)) AS hemoglobin_high_count,
  MAX(IF(signal = 'hemoglobin', bin_low_count, NULL)) AS hemoglobin_bin_low,
  MAX(IF(signal = 'hemoglobin', bin_normal_count, NULL)) AS hemoglobin_bin_normal,
  MAX(IF(signal = 'hemoglobin', bin_high_count, NULL)) AS hemoglobin_bin_high,

  -- =========================
  -- GLUCOSE
  -- =========================
  MAX(IF(signal = 'glucose', value_min, NULL)) AS glucose_min,
  MAX(IF(signal = 'glucose', value_max, NULL)) AS glucose_max,
  MAX(IF(signal = 'glucose', value_mean, NULL)) AS glucose_mean,
  MAX(IF(signal = 'glucose', value_median, NULL)) AS glucose_median,
  MAX(IF(signal = 'glucose', value_iqr, NULL)) AS glucose_iqr,
  MAX(IF(signal = 'glucose', measurement_count, NULL)) AS glucose_count,
  MAX(IF(signal = 'glucose', missing_indicator, NULL)) AS glucose_missing,
  MAX(IF(signal = 'glucose', early_mean, NULL)) AS glucose_early_mean,
  MAX(IF(signal = 'glucose', late_mean, NULL)) AS glucose_late_mean,
  MAX(IF(signal = 'glucose', count_below_safe_min, NULL)) AS glucose_low_count,
  MAX(IF(signal = 'glucose', count_above_safe_max, NULL)) AS glucose_high_count,
  MAX(IF(signal = 'glucose', bin_low_count, NULL)) AS glucose_bin_low,
  MAX(IF(signal = 'glucose', bin_normal_count, NULL)) AS glucose_bin_normal,
  MAX(IF(signal = 'glucose', bin_high_count, NULL)) AS glucose_bin_high,

  -- =========================
  -- WBC
  -- =========================
  MAX(IF(signal = 'wbc', value_min, NULL)) AS wbc_min,
  MAX(IF(signal = 'wbc', value_max, NULL)) AS wbc_max,
  MAX(IF(signal = 'wbc', value_mean, NULL)) AS wbc_mean,
  MAX(IF(signal = 'wbc', value_median, NULL)) AS wbc_median,
  MAX(IF(signal = 'wbc', value_iqr, NULL)) AS wbc_iqr,
  MAX(IF(signal = 'wbc', measurement_count, NULL)) AS wbc_count,
  MAX(IF(signal = 'wbc', missing_indicator, NULL)) AS wbc_missing,
  MAX(IF(signal = 'wbc', early_mean, NULL)) AS wbc_early_mean,
  MAX(IF(signal = 'wbc', late_mean, NULL)) AS wbc_late_mean,
  MAX(IF(signal = 'wbc', count_below_safe_min, NULL)) AS wbc_low_count,
  MAX(IF(signal = 'wbc', count_above_safe_max, NULL)) AS wbc_high_count,
  MAX(IF(signal = 'wbc', bin_low_count, NULL)) AS wbc_bin_low,
  MAX(IF(signal = 'wbc', bin_normal_count, NULL)) AS wbc_bin_normal,
  MAX(IF(signal = 'wbc', bin_high_count, NULL)) AS wbc_bin_high,

  -- =========================
  -- PLATELET
  -- =========================
  MAX(IF(signal = 'platelet', value_min, NULL)) AS platelet_min,
  MAX(IF(signal = 'platelet', value_max, NULL)) AS platelet_max,
  MAX(IF(signal = 'platelet', value_mean, NULL)) AS platelet_mean,
  MAX(IF(signal = 'platelet', value_median, NULL)) AS platelet_median,
  MAX(IF(signal = 'platelet', value_iqr, NULL)) AS platelet_iqr,
  MAX(IF(signal = 'platelet', measurement_count, NULL)) AS platelet_count,
  MAX(IF(signal = 'platelet', missing_indicator, NULL)) AS platelet_missing,
  MAX(IF(signal = 'platelet', early_mean, NULL)) AS platelet_early_mean,
  MAX(IF(signal = 'platelet', late_mean, NULL)) AS platelet_late_mean,
  MAX(IF(signal = 'platelet', count_below_safe_min, NULL)) AS platelet_low_count,
  MAX(IF(signal = 'platelet', count_above_safe_max, NULL)) AS platelet_high_count,
  MAX(IF(signal = 'platelet', bin_low_count, NULL)) AS platelet_bin_low,
  MAX(IF(signal = 'platelet', bin_normal_count, NULL)) AS platelet_bin_normal,
  MAX(IF(signal = 'platelet', bin_high_count, NULL)) AS platelet_bin_high,

  -- =========================
  -- LACTATE
  -- =========================
  MAX(IF(signal = 'lactate', value_min, NULL)) AS lactate_min,
  MAX(IF(signal = 'lactate', value_max, NULL)) AS lactate_max,
  MAX(IF(signal = 'lactate', value_mean, NULL)) AS lactate_mean,
  MAX(IF(signal = 'lactate', value_median, NULL)) AS lactate_median,
  MAX(IF(signal = 'lactate', value_iqr, NULL)) AS lactate_iqr,
  MAX(IF(signal = 'lactate', measurement_count, NULL)) AS lactate_count,
  MAX(IF(signal = 'lactate', missing_indicator, NULL)) AS lactate_missing,
  MAX(IF(signal = 'lactate', early_mean, NULL)) AS lactate_early_mean,
  MAX(IF(signal = 'lactate', late_mean, NULL)) AS lactate_late_mean,
  MAX(IF(signal = 'lactate', count_below_safe_min, NULL)) AS lactate_low_count,
  MAX(IF(signal = 'lactate', count_above_safe_max, NULL)) AS lactate_high_count,
  MAX(IF(signal = 'lactate', bin_low_count, NULL)) AS lactate_bin_low,
  MAX(IF(signal = 'lactate', bin_normal_count, NULL)) AS lactate_bin_normal,
  MAX(IF(signal = 'lactate', bin_high_count, NULL)) AS lactate_bin_high

FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_24h_signal_features_long_2`
GROUP BY stay_id, mortality, gender, anchor_age;