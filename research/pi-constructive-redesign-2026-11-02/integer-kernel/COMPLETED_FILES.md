# Completed evidence

Main theorem and constructive direction:
- INTEGER_ANTIDERIVATIVE_FAMILY.md
- verify_family.py / family_checks.json
- verify_odd_period_image.py / fair_endpoint_comparison.json

Finite construction/search controls:
- antiderivative_search.py / antiderivative_search_results.json (initial integer sublattice, through c=20)
- antiderivative_small_saturated.py / small_saturated_results.json (full integer kernel, through c=12)
- compare_single_kernels.py / single_kernel_control.json
- verify_constructive_bounds.py / constructive_upper_certificates.json
- coordinate_cost_table.json and period_lattice_controls.json (descriptive controls)

The large antiderivative_search_saturated.py c<=32 run was stopped and has no completed full-range result. It remains source code for future bounded work, not completed evidence. No claim of a new rational approximation record or an infinite quality-zero sequence is made.

Additional completed fixed growing family:
- new_family_search.py / new_family_results.json: m=1,...,6 all completed in about 8.53 seconds, a=m,b=2m,c=4m+1,H_degree<=m+1; all full coordinates and exact quality intervals included.
