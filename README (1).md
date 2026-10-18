# Fixed primary-model equation

Coefficients and reference survival were exported from the final 100-mm clinical-only Cox model (1,243 patients, 310 deaths). Each estimate is based on the regression's reference design values.

\`LP = Σ beta_j (x_j − x_ref,j)\`

\`P(death by t) = 1 − S_ref(t)^exp(LP)\`, where \`t\` is 60 or 120 months.

- \`AGE_NUM\` is age in years.
- \`tumor_size_10mm\` is maximal recorded tumor size in cm (equivalently millimeters divided by 10).
- \`SEX3Female\` equals 1 for female and 0 for male.
- \`CDCC21\` equals 1 for Charlson–Deyo category 1.
- \`CDCC22+\` equals 1 for Charlson–Deyo category ≥2.
- Charlson–Deyo category 0 is the reference group.

The coefficients are *fixed manuscript outputs*, not sample data. Do not replace the validated equation with coefficients fitted on an unrelated cohort. Internal validation does not establish suitability for clinical decision-making.
