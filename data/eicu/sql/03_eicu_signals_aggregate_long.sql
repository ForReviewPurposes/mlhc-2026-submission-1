CREATE OR REPLACE TABLE `eicu-ext.eicu_ext_data.eicu_signal_24h_agg_long` AS
WITH vals AS (
  SELECT * FROM `eicu-ext.eicu_ext_data.eicu_vals_raw`
),
cohort AS (
  SELECT * FROM `eicu-ext.eicu_ext_data.eicu_cohort_all`
),
signal_ranges AS (
  SELECT 'heart_rate' AS signal, 50.0 AS safe_min, 110.0 AS safe_max
  UNION ALL SELECT 'sbp', 90.0, 180.0
  UNION ALL SELECT 'dbp', 50.0, 100.0
  UNION ALL SELECT 'mbp', 65.0, 110.0
  UNION ALL SELECT 'resp_rate', 10.0, 25.0
  UNION ALL SELECT 'spo2', 92.0, 100.0
  UNION ALL SELECT 'temperature', 36.0, 38.3
  UNION ALL SELECT 'creatinine', 0.6, 1.3
  UNION ALL SELECT 'glucose', 70.0, 180.0
  UNION ALL SELECT 'sodium', 135.0, 145.0
  UNION ALL SELECT 'potassium', 3.5, 5.2
  UNION ALL SELECT 'bicarbonate', 22.0, 30.0
  UNION ALL SELECT 'wbc', 4.0, 11.0
  UNION ALL SELECT 'hemoglobin', 10.0, 17.0
  UNION ALL SELECT 'platelet', 150.0, 400.0
  UNION ALL SELECT 'lactate', 0.5, 2.0
  UNION ALL SELECT 'gcs', 13.0, 15.0
  UNION ALL SELECT 'gcs_motor', 5.0, 6.0
  UNION ALL SELECT 'gcs_verbal', 4.0, 5.0
  UNION ALL SELECT 'gcs_eyes', 3.0, 4.0
),
stay_signal_grid AS (
  SELECT
    c.patientunitstayid,
    c.hospitalid,
    sr.signal,
    sr.safe_min,
    sr.safe_max
  FROM cohort c
  CROSS JOIN signal_ranges sr
),
agg AS (
  SELECT
    g.patientunitstayid,
    g.hospitalid,
    c.mortality,
    c.gender,
    c.anchor_age,
    c.duration_hours_from_24h,
    c.event_after_24h,
    g.signal,
    g.safe_min,
    g.safe_max,
    COUNT(v.value) AS measurement_count,
    IF(COUNT(v.value) = 0, 1, 0) AS missing_indicator,
    MIN(v.value) AS value_min,
    MAX(v.value) AS value_max,
    AVG(v.value) AS value_mean,
    APPROX_QUANTILES(v.value, 4)[OFFSET(1)] AS value_q25,
    APPROX_QUANTILES(v.value, 4)[OFFSET(2)] AS value_median,
    APPROX_QUANTILES(v.value, 4)[OFFSET(3)] AS value_q75,
    AVG(IF(v.offset < 720, v.value, NULL)) AS early_mean,
    AVG(IF(v.offset >= 720, v.value, NULL)) AS late_mean,
    COUNTIF(v.value < g.safe_min) AS count_below_safe_min,
    COUNTIF(v.value > g.safe_max) AS count_above_safe_max,
    COUNTIF(v.value < g.safe_min) AS bin_low_count,
    COUNTIF(v.value >= g.safe_min AND v.value <= g.safe_max) AS bin_normal_count,
    COUNTIF(v.value > g.safe_max) AS bin_high_count
  FROM stay_signal_grid g
  JOIN cohort c
    ON g.patientunitstayid = c.patientunitstayid
  LEFT JOIN vals v
    ON g.patientunitstayid = v.patientunitstayid
   AND g.signal = v.signal
  GROUP BY
    g.patientunitstayid,
    g.hospitalid,
    c.mortality,
    c.gender,
    c.anchor_age,
    c.duration_hours_from_24h,
    c.event_after_24h,
    g.signal,
    g.safe_min,
    g.safe_max
)
SELECT
  *,
  value_q75 - value_q25 AS value_iqr
FROM agg;