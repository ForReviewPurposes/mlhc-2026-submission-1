-- Query 1: GCS long format extraction
-- Follows the same structure as mimiciv_ext_icumort_signal24h_extract_1.sql
-- GCS total: 3-15, normal = 15 (fully alert)
-- GCS motor: 1-6, normal = 6
-- GCS verbal: 1-5, normal = 5
-- GCS eyes: 1-4, normal = 4
-- Safe ranges: "normal" is max score; lower scores indicate impairment
-- For GCS, safe_min is set at thresholds below which concern increases

CREATE OR REPLACE TABLE `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_24h_gcs_features_long_2`
AS
WITH
  signal_ranges AS (
    SELECT 'gcs' AS signal, 13.0 AS safe_min, 15.0 AS safe_max
    UNION ALL
    SELECT 'gcs_motor', 5.0, 6.0
    UNION ALL
    SELECT 'gcs_verbal', 4.0, 5.0
    UNION ALL
    SELECT 'gcs_eyes', 3.0, 4.0
  ),
  cohort AS (
    SELECT *
    FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_icu_cohort_1_timetodeath`
  ),
  vals AS (
    -- GCS total
    SELECT
      c.stay_id,
      'gcs' AS signal,
      g.charttime,
      CAST(g.gcs AS FLOAT64) AS value
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.gcs` g
      ON
        c.stay_id = g.stay_id
        AND g.charttime >= c.intime
        AND g.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE g.gcs IS NOT NULL

    UNION ALL

    -- GCS motor
    SELECT
      c.stay_id,
      'gcs_motor' AS signal,
      g.charttime,
      CAST(g.gcs_motor AS FLOAT64) AS value
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.gcs` g
      ON
        c.stay_id = g.stay_id
        AND g.charttime >= c.intime
        AND g.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE g.gcs_motor IS NOT NULL

    UNION ALL

    -- GCS verbal
    SELECT
      c.stay_id,
      'gcs_verbal' AS signal,
      g.charttime,
      CAST(g.gcs_verbal AS FLOAT64) AS value
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.gcs` g
      ON
        c.stay_id = g.stay_id
        AND g.charttime >= c.intime
        AND g.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE g.gcs_verbal IS NOT NULL

    UNION ALL

    -- GCS eyes
    SELECT
      c.stay_id,
      'gcs_eyes' AS signal,
      g.charttime,
      CAST(g.gcs_eyes AS FLOAT64) AS value
    FROM cohort c
    JOIN `physionet-data.mimiciv_3_1_derived.gcs` g
      ON
        c.stay_id = g.stay_id
        AND g.charttime >= c.intime
        AND g.charttime < TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
    WHERE g.gcs_eyes IS NOT NULL
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
