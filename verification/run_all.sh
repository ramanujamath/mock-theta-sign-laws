#!/bin/sh
# Regenerates the data files and runs every verification script, saving outputs in results/.
# Usage: sh run_all.sh        (several hours on 2 cores; individual scripts can be run directly)
cd "$(dirname "$0")"
mkdir -p results
run() { name=$1; shift; echo "== $name $*"; start=$(date +%s)
        timeout 3000 python3 "$@" > "results/$name.txt" 2>&1; rc=$?
        echo "   exit $rc, $(( $(date +%s) - start )) s"; }
# data generation (pickles, not committed)
run big          big.py 60000        # exact coefficients S0, S1, mu2 to n = 60000 -> big.pkl
run family       family.py 20000     # exact coefficients U0, U1, T0, T1, V0, V1 to n = 20000 -> family.pkl
run mock         mock.py             # third-order phi, chi to n = 1500 -> mt.pkl
run scan         scan.py             # other orders, sign-period scan -> F.pkl
run th_check     th_check.py         # TH coefficients -> TH.pkl; identities U0+2U1=TH, mu2=2U0-TH to q^20000
# finite checks and identities
run review_finite review_finite.py
run review_series review_series.py
run verify       verify.py
run repcheck2    repcheck2.py
run al           al.py
run altforms     altforms.py
run review_forms review_forms.py
run prod_identity prod_identity.py
run exact        exact.py
run exact_phase  exact_phase.py
run th_const     th_const.py
run review_phase review_phase.py
# asymptotics and radial behaviour
run radial       radial.py
run pref         pref.py
run fit3         fit3.py
run scan2        scan2.py
run split        split.py
run resid        resid.py
run u0           u0.py
run u0b          u0b.py
run ufit         ufit.py
# cusp exponents
run cusps        cusps.py
run cusps2       cusps2.py
run lemma        lemma.py
run farey_exponents farey_exponents.py
run farey_S1_exp farey_S1_exp.py
run review_exponents review_exponents.py
run review_cmpR  review_cmpR.py
run u_cusps      u_cusps.py search exact cusps radial
# circle-method bounds and their validation
run farey_check_transform farey_check_transform.py
run farey_validate farey_validate.py
run farey_S1_validate farey_S1_validate.py
run farey_S1_v1dbg farey_S1_v1dbg.py
run farey_gen_validate farey_gen_validate.py
run farey_gen_v1b farey_gen_v1b.py
run farey_gen_v5 farey_gen_v5.py
run review_circle review_circle.py
run review_circle2 review_circle2.py
run review_center review_center.py
run review_major_S0 review_major.py S0 3000
run review_major_S1 review_major.py S1 3000
run review_verify review_verify.py
run review_cert_S0 review_cert.py S0 11200 11300
run review_cert_mu2 review_cert.py mu2 950 1100
run certify_run  certify_run.py
