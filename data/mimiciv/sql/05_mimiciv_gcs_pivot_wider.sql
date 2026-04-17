-- Query 2: GCS wide format pivot
-- Follows the same structure as mimiciv_ext_icumort_24h_features_long_extract_2.sql
-- Produces one row per stay_id with GCS columns to join with existing wide table

CREATE OR REPLACE TABLE `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_24h_gcs_features_wide_2`
AS
SELECT
  stay_id,

  -- =========================
  -- GCS TOTAL
  -- =========================
  MAX(IF(signal = 'gcs', value_min, NULL)) AS gcs_min,
  MAX(IF(signal = 'gcs', value_max, NULL)) AS gcs_max,
  MAX(IF(signal = 'gcs', value_mean, NULL)) AS gcs_mean,
  MAX(IF(signal = 'gcs', value_median, NULL)) AS gcs_median,
  MAX(IF(signal = 'gcs', value_iqr, NULL)) AS gcs_iqr,
  MAX(IF(signal = 'gcs', measurement_count, NULL)) AS gcs_count,
  MAX(IF(signal = 'gcs', missing_indicator, NULL)) AS gcs_missing,
  MAX(IF(signal = 'gcs', early_mean, NULL)) AS gcs_early_mean,
  MAX(IF(signal = 'gcs', late_mean, NULL)) AS gcs_late_mean,
  MAX(IF(signal = 'gcs', count_below_safe_min, NULL)) AS gcs_low_count,
  MAX(IF(signal = 'gcs', count_above_safe_max, NULL)) AS gcs_high_count,
  MAX(IF(signal = 'gcs', bin_low_count, NULL)) AS gcs_bin_low,
  MAX(IF(signal = 'gcs', bin_normal_count, NULL)) AS gcs_bin_normal,
  MAX(IF(signal = 'gcs', bin_high_count, NULL)) AS gcs_bin_high,

  -- =========================
  -- GCS MOTOR
  -- =========================
  MAX(IF(signal = 'gcs_motor', value_min, NULL)) AS gcs_motor_min,
  MAX(IF(signal = 'gcs_motor', value_max, NULL)) AS gcs_motor_max,
  MAX(IF(signal = 'gcs_motor', value_mean, NULL)) AS gcs_motor_mean,
  MAX(IF(signal = 'gcs_motor', value_median, NULL)) AS gcs_motor_median,
  MAX(IF(signal = 'gcs_motor', value_iqr, NULL)) AS gcs_motor_iqr,
  MAX(IF(signal = 'gcs_motor', measurement_count, NULL)) AS gcs_motor_count,
  MAX(IF(signal = 'gcs_motor', missing_indicator, NULL)) AS gcs_motor_missing,
  MAX(IF(signal = 'gcs_motor', early_mean, NULL)) AS gcs_motor_early_mean,
  MAX(IF(signal = 'gcs_motor', late_mean, NULL)) AS gcs_motor_late_mean,
  MAX(IF(signal = 'gcs_motor', count_below_safe_min, NULL)) AS gcs_motor_low_count,
  MAX(IF(signal = 'gcs_motor', count_above_safe_max, NULL)) AS gcs_motor_high_count,
  MAX(IF(signal = 'gcs_motor', bin_low_count, NULL)) AS gcs_motor_bin_low,
  MAX(IF(signal = 'gcs_motor', bin_normal_count, NULL)) AS gcs_motor_bin_normal,
  MAX(IF(signal = 'gcs_motor', bin_high_count, NULL)) AS gcs_motor_bin_high,

  -- =========================
  -- GCS VERBAL
  -- =========================
  MAX(IF(signal = 'gcs_verbal', value_min, NULL)) AS gcs_verbal_min,
  MAX(IF(signal = 'gcs_verbal', value_max, NULL)) AS gcs_verbal_max,
  MAX(IF(signal = 'gcs_verbal', value_mean, NULL)) AS gcs_verbal_mean,
  MAX(IF(signal = 'gcs_verbal', value_median, NULL)) AS gcs_verbal_median,
  MAX(IF(signal = 'gcs_verbal', value_iqr, NULL)) AS gcs_verbal_iqr,
  MAX(IF(signal = 'gcs_verbal', measurement_count, NULL)) AS gcs_verbal_count,
  MAX(IF(signal = 'gcs_verbal', missing_indicator, NULL)) AS gcs_verbal_missing,
  MAX(IF(signal = 'gcs_verbal', early_mean, NULL)) AS gcs_verbal_early_mean,
  MAX(IF(signal = 'gcs_verbal', late_mean, NULL)) AS gcs_verbal_late_mean,
  MAX(IF(signal = 'gcs_verbal', count_below_safe_min, NULL)) AS gcs_verbal_low_count,
  MAX(IF(signal = 'gcs_verbal', count_above_safe_max, NULL)) AS gcs_verbal_high_count,
  MAX(IF(signal = 'gcs_verbal', bin_low_count, NULL)) AS gcs_verbal_bin_low,
  MAX(IF(signal = 'gcs_verbal', bin_normal_count, NULL)) AS gcs_verbal_bin_normal,
  MAX(IF(signal = 'gcs_verbal', bin_high_count, NULL)) AS gcs_verbal_bin_high,

  -- =========================
  -- GCS EYES
  -- =========================
  MAX(IF(signal = 'gcs_eyes', value_min, NULL)) AS gcs_eyes_min,
  MAX(IF(signal = 'gcs_eyes', value_max, NULL)) AS gcs_eyes_max,
  MAX(IF(signal = 'gcs_eyes', value_mean, NULL)) AS gcs_eyes_mean,
  MAX(IF(signal = 'gcs_eyes', value_median, NULL)) AS gcs_eyes_median,
  MAX(IF(signal = 'gcs_eyes', value_iqr, NULL)) AS gcs_eyes_iqr,
  MAX(IF(signal = 'gcs_eyes', measurement_count, NULL)) AS gcs_eyes_count,
  MAX(IF(signal = 'gcs_eyes', missing_indicator, NULL)) AS gcs_eyes_missing,
  MAX(IF(signal = 'gcs_eyes', early_mean, NULL)) AS gcs_eyes_early_mean,
  MAX(IF(signal = 'gcs_eyes', late_mean, NULL)) AS gcs_eyes_late_mean,
  MAX(IF(signal = 'gcs_eyes', count_below_safe_min, NULL)) AS gcs_eyes_low_count,
  MAX(IF(signal = 'gcs_eyes', count_above_safe_max, NULL)) AS gcs_eyes_high_count,
  MAX(IF(signal = 'gcs_eyes', bin_low_count, NULL)) AS gcs_eyes_bin_low,
  MAX(IF(signal = 'gcs_eyes', bin_normal_count, NULL)) AS gcs_eyes_bin_normal,
  MAX(IF(signal = 'gcs_eyes', bin_high_count, NULL)) AS gcs_eyes_bin_high

FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_24h_gcs_features_long_2`
GROUP BY stay_id;
