CREATE OR REPLACE TABLE `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_24h_signal_features_long_2`
AS
WITH
  signal_ranges AS (
    SELECT 'heart_rate' AS signal, 50.0 AS safe_min, 110.0 AS safe_max
    UNION ALL
    SELECT 'sbp', 90.0, 180.0
    UNION ALL
    SELECT 'dbp', 50.0, 100.0
    UNION ALL
    SELECT 'mbp', 65.0, 110.0
    UNION ALL
    SELECT 'resp_rate', 10.0, 25.0
    UNION ALL
    SELECT 'spo2', 92.0, 100.0
    UNION ALL
    SELECT 'temperature', 36.0, 38.3
    UNION ALL
    SELECT 'creatinine', 0.6, 1.3
    UNION ALL
    SELECT 'glucose', 70.0, 180.0
    UNION ALL
    SELECT 'sodium', 135.0, 145.0
    UNION ALL
    SELECT 'potassium', 3.5, 5.2
    UNION ALL
    SELECT 'bicarbonate', 22.0, 30.0
    UNION ALL
    SELECT 'wbc', 4.0, 11.0
    UNION ALL
    SELECT 'hemoglobin', 10.0, 17.0
    UNION ALL
    SELECT 'platelet', 150.0, 400.0
    UNION ALL
    SELECT 'lactate', 0.5, 2.0
  ),
  cohort AS (
    SELECT *
    FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_icu_cohort_1_timetodeath`
  ),
  vals AS (
    -- vitals
    SELECT
      c.stay_id,
      'heart_rate' AS signal,
      v.charttime,
      CAST(v.heart_rate AS FLOAT64) AS value
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
      ON
        c.stay_id = v.stay_id
        AND v.charttime >= c.intime
        AND v.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE v.heart_rate IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'sbp', v.charttime, CAST(v.sbp AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
      ON
        c.stay_id = v.stay_id
        AND v.charttime >= c.intime
        AND v.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE v.sbp IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'dbp', v.charttime, CAST(v.dbp AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
      ON
        c.stay_id = v.stay_id
        AND v.charttime >= c.intime
        AND v.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE v.dbp IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'mbp', v.charttime, CAST(v.mbp AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
      ON
        c.stay_id = v.stay_id
        AND v.charttime >= c.intime
        AND v.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE v.mbp IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'resp_rate', v.charttime, CAST(v.resp_rate AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
      ON
        c.stay_id = v.stay_id
        AND v.charttime >= c.intime
        AND v.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE v.resp_rate IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'spo2', v.charttime, CAST(v.spo2 AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
      ON
        c.stay_id = v.stay_id
        AND v.charttime >= c.intime
        AND v.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE v.spo2 IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'temperature', v.charttime, CAST(v.temperature AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
      ON
        c.stay_id = v.stay_id
        AND v.charttime >= c.intime
        AND v.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE v.temperature IS NOT NULL

    -- chemistry
    UNION ALL
    SELECT c.stay_id, 'creatinine', ch.charttime, CAST(ch.creatinine AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
      ON
        c.hadm_id = ch.hadm_id
        AND ch.charttime >= c.intime
        AND ch.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE ch.creatinine IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'glucose', ch.charttime, CAST(ch.glucose AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
      ON
        c.hadm_id = ch.hadm_id
        AND ch.charttime >= c.intime
        AND ch.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE ch.glucose IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'sodium', ch.charttime, CAST(ch.sodium AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
      ON
        c.hadm_id = ch.hadm_id
        AND ch.charttime >= c.intime
        AND ch.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE ch.sodium IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'potassium', ch.charttime, CAST(ch.potassium AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
      ON
        c.hadm_id = ch.hadm_id
        AND ch.charttime >= c.intime
        AND ch.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE ch.potassium IS NOT NULL
    UNION ALL
    SELECT
      c.stay_id, 'bicarbonate', ch.charttime, CAST(ch.bicarbonate AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
      ON
        c.hadm_id = ch.hadm_id
        AND ch.charttime >= c.intime
        AND ch.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE ch.bicarbonate IS NOT NULL

    -- CBC
    UNION ALL
    SELECT c.stay_id, 'wbc', cbc.charttime, CAST(cbc.wbc AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.complete_blood_count` cbc
      ON
        c.hadm_id = cbc.hadm_id
        AND cbc.charttime >= c.intime
        AND cbc.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE cbc.wbc IS NOT NULL
    UNION ALL
    SELECT
      c.stay_id, 'hemoglobin', cbc.charttime, CAST(cbc.hemoglobin AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.complete_blood_count` cbc
      ON
        c.hadm_id = cbc.hadm_id
        AND cbc.charttime >= c.intime
        AND cbc.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE cbc.hemoglobin IS NOT NULL
    UNION ALL
    SELECT c.stay_id, 'platelet', cbc.charttime, CAST(cbc.platelet AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.complete_blood_count` cbc
      ON
        c.hadm_id = cbc.hadm_id
        AND cbc.charttime >= c.intime
        AND cbc.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE cbc.platelet IS NOT NULL

    -- blood gas
    UNION ALL
    SELECT c.stay_id, 'lactate', bg.charttime, CAST(bg.lactate AS FLOAT64)
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.bg` bg
      ON
        c.hadm_id = bg.hadm_id
        AND bg.charttime >= c.intime
        AND bg.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE bg.lactate IS NOT NULL
  ),
  stay_signal_grid AS (
    SELECT
      c.stay_id,
      c.mortality,
      c.gender,
      c.anchor_age,
      sr.signal,
      sr.safe_min,
      sr.safe_max
    FROM cohort c
    CROSS JOIN signal_ranges sr
  ),
  agg AS (
    SELECT
      g.stay_id,
      g.mortality,
      g.gender,
      g.anchor_age,
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
      AVG(
        IF(
          v.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 12 HOUR),
          v.value,
          NULL)) AS early_mean,
      AVG(
        IF(
          v.charttime >= TIMESTAMP_ADD(c.intime, INTERVAL 12 HOUR),
          v.value,
          NULL)) AS late_mean,
      COUNTIF(v.value < g.safe_min) AS count_below_safe_min,
      COUNTIF(v.value > g.safe_max) AS count_above_safe_max,
      COUNTIF(v.value < g.safe_min) AS bin_low_count,
      COUNTIF(v.value >= g.safe_min AND v.value <= g.safe_max)
        AS bin_normal_count,
      COUNTIF(v.value > g.safe_max) AS bin_high_count
    FROM stay_signal_grid g
    JOIN cohort c
      ON g.stay_id = c.stay_id
    LEFT JOIN vals v
      ON
        g.stay_id = v.stay_id
        AND g.signal = v.signal
    GROUP BY
      g.stay_id, g.mortality, g.gender, g.anchor_age, g.signal, g.safe_min,
      g.safe_max, c.intime
  )
SELECT
  *,
  value_q75 - value_q25 AS value_iqr
FROM agg;
