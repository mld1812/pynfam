# Default settings for one fit
# They overwrite the default settings given in PyNFAM
# They can be overwritten by "setts_fit" in pynfam_fit
default_setts = {
    'pynfam_inputs':  {
        'gs_def_scan': (-2, (-0.2,0.0,0.2)),
        'beta_type': '-',
        'fam_ops': 'ALL'
    },
    'hfb':  {
        'number_of_shells': 20,
        'accuracy': 1e-06,
        'restart_file': 1,
        'number_iterations': 300,
        'functional': 'HFB1',
        'user_pairing': False,
        'vpair_n': 0.0,
        'vpair_p': 0.0,
        'pairing_cutoff': 60.0,
        'pairing_feature': 0.5
            },
    'fam':  {
        'use_fam_storage': 0,
        'interaction_name': 'HFB1',
        'force_j2_terms': False,
        'convergence_epsilon': 1e-07,
        'require_self_consistency': True,
        'require_gauge_invariance': False,
        'vpair_t0': 0.0,
        'vpair_t1': None,
        'override_cs0': None,
        'override_csr': 0.0,
        'override_cds': 0.0,
        'override_ct': 0.0,
        'override_cgs': 0.0,
        'override_cf': 0.0,
        'override_cj': None,
        'override_csdj': None,
        'max_iter': 300,
        'broyden_history_size': 170
            },
    'ctr':  {
        'half_width': 0.1, # for resonance calc only
        'de_hw_ratio': 1.0, # for resonance calc only
        'nr_points': 60 # for half-live calc only
            },
    'psi':  {
        'GV': 1.0,
        'GA': -1.0
            },
    'scaled_params': {
        'x1': 1.0, # Cs has no density dependence
        'g0p': 1.6, # Sn132 experiment result
        'g1p': 0.0, # no tensor term
        'h0p': 0.0, # no tensor term
        # isoscalar pairing scaled w.r.t. the absolute value of the average isovector pairing
        'vpair_t0_scaled': 0.0
            }
    }