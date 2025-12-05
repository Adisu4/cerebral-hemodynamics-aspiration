# Imports moved to top for best practice
import numpy as np
import matplotlib.pyplot as plt
def stenosis_patterns_pvs_analysis():
	"""
	Analyze the effects of different venous stenosis patterns (A, B, C, D) on intracranial pressure (ICP) and venous sinus pressure (Pvs).
	Simulates both complete blockage (Gx=0) and halved conductance (Gx/2) for each pattern in supine and upright postures.
	Generates bar plots and summary tables for all results.
	"""
	print("\n" + "="*70)
	print("STENOTIC PATTERNS ANALYSIS: ICP and Pvs")
	print("="*70)
	patterns = ['Baseline', 'Pattern A', 'Pattern B', 'Pattern C', 'Pattern D']
	icp_supine_blocked = []
	icp_upright_blocked = []
	icp_supine_halved = []
	icp_upright_halved = []
	pvs_supine_blocked = []
	pvs_upright_blocked = []
	pvs_supine_halved = []
	pvs_upright_halved = []
	y0 = ModelParameters().get_initial_conditions()
	for pattern in patterns:
		print(f"\n{'='*70}")
		print(f"Pattern: {pattern}")
		print(f"{'='*70}")
		# COMPLETE BLOCKAGE (Gx = 0)
		print(f"\n{pattern} - COMPLETE BLOCKAGE (Gx=0)")
		print("-"*70)
		params_sup_block = ModelParameters(posture='supine')
		params_upr_block = ModelParameters(posture='upright')
		# Apply pattern
		if pattern == 'Pattern A':
			params_sup_block.Gazy2 = 0.0
			params_sup_block.kjl3 = 0.0
			params_sup_block.kjl2 = 0.0
			params_sup_block.kjl1 = 0.0
			params_upr_block.Gazy2 = 0.0
			params_upr_block.kjl3 = 0.0
			params_upr_block.kjl2 = 0.0
			params_upr_block.kjl1 = 0.0
		elif pattern == 'Pattern B':
			for p in [params_sup_block, params_upr_block]:
				p.kjr3 = 0.0; p.kjl3 = 0.0; p.kjr2 = 0.0; p.kjl2 = 0.0; p.kjr1 = 0.0; p.kjl1 = 0.0; p.Gazy2 = 0.0
		elif pattern == 'Pattern C':
			for p in [params_sup_block, params_upr_block]:
				p.kjr3 = 0.0; p.kjl3 = 0.0; p.kjr2 = 0.0; p.kjl2 = 0.0; p.kjr1 = 0.0; p.kjl1 = 0.0
		elif pattern == 'Pattern D':
			for p in [params_sup_block, params_upr_block]:
				p.Gazy1 = 0.0; p.Gazy2 = 0.0; p.Glv = 0.0
		# SUPINE
		sol = solve_ivp(gadda_ode_system_FIXED, [0, 500], y0, args=(params_sup_block,),
					   method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8)
		if sol.success:
			icp_supine_blocked.append(sol.y[0, -1])
			pvs_supine_blocked.append(sol.y[3, -1])
		else:
			icp_supine_blocked.append(np.nan)
			pvs_supine_blocked.append(np.nan)
		# UPRIGHT
		sol = solve_ivp(gadda_ode_system_FIXED, [0, 500], y0, args=(params_upr_block,),
					   method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8)
		if sol.success:
			icp_upright_blocked.append(sol.y[0, -1])
			pvs_upright_blocked.append(sol.y[3, -1])
		else:
			icp_upright_blocked.append(np.nan)
			pvs_upright_blocked.append(np.nan)
		# HALVED CONDUCTANCE (Gx / 2)
		print(f"\n{pattern} - HALVED CONDUCTANCE (Gx/2)")
		print("-"*70)
		params_sup_half = ModelParameters(posture='supine')
		params_upr_half = ModelParameters(posture='upright')
		if pattern == 'Pattern A':
			params_sup_half.Gazy2 /= 2.0; params_sup_half.kjl3 /= 2.0; params_sup_half.kjl2 /= 2.0; params_sup_half.kjl1 /= 2.0
			params_upr_half.Gazy2 /= 2.0; params_upr_half.kjl3 /= 2.0; params_upr_half.kjl2 /= 2.0; params_upr_half.kjl1 /= 2.0
		elif pattern == 'Pattern B':
			for p in [params_sup_half, params_upr_half]:
				p.kjr3 /= 2.0; p.kjl3 /= 2.0; p.kjr2 /= 2.0; p.kjl2 /= 2.0; p.kjr1 /= 2.0; p.kjl1 /= 2.0; p.Gazy2 /= 2.0
		elif pattern == 'Pattern C':
			for p in [params_sup_half, params_upr_half]:
				p.kjr3 /= 2.0; p.kjl3 /= 2.0; p.kjr2 /= 2.0; p.kjl2 /= 2.0; p.kjr1 /= 2.0; p.kjl1 /= 2.0
		elif pattern == 'Pattern D':
			for p in [params_sup_half, params_upr_half]:
				p.Gazy1 /= 2.0; p.Gazy2 /= 2.0; p.Glv /= 2.0
		# SUPINE
		sol = solve_ivp(gadda_ode_system_FIXED, [0, 500], y0, args=(params_sup_half,),
					   method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8)
		if sol.success:
			icp_supine_halved.append(sol.y[0, -1])
			pvs_supine_halved.append(sol.y[3, -1])
		else:
			icp_supine_halved.append(np.nan)
			pvs_supine_halved.append(np.nan)
		# UPRIGHT
		sol = solve_ivp(gadda_ode_system_FIXED, [0, 500], y0, args=(params_upr_half,),
					   method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8)
		if sol.success:
			icp_upright_halved.append(sol.y[0, -1])
			pvs_upright_halved.append(sol.y[3, -1])
		else:
			icp_upright_halved.append(np.nan)
			pvs_upright_halved.append(np.nan)
	# Print summary
	print("\n" + "="*70)
	print("STENOTIC PATTERNS SUMMARY - ICP (mmHg)")
	print("="*70)
	print(f"{'Pattern':<15} {'Sup Gx=0':<12} {'Upr Gx=0':<12} {'Sup Gx/2':<12} {'Upr Gx/2':<12}")
	print("-"*70)
	for i, pattern in enumerate(patterns):
		print(f"{pattern:<15} {icp_supine_blocked[i]:>10.3f}  {icp_upright_blocked[i]:>10.3f}  "
			  f"{icp_supine_halved[i]:>10.3f}  {icp_upright_halved[i]:>10.3f}")
	print("="*70)
	print("\n" + "="*70)
	print("STENOTIC PATTERNS SUMMARY - Pvs (mmHg)")
	print("="*70)
	print(f"{'Pattern':<15} {'Sup Gx=0':<12} {'Upr Gx=0':<12} {'Sup Gx/2':<12} {'Upr Gx/2':<12}")
	print("-"*70)
	for i, pattern in enumerate(patterns):
		print(f"{pattern:<15} {pvs_supine_blocked[i]:>10.3f}  {pvs_upright_blocked[i]:>10.3f}  "
			  f"{pvs_supine_halved[i]:>10.3f}  {pvs_upright_halved[i]:>10.3f}")
	print("="*70)
	# Plotting: Bar plots for ICP and Pvs (supine/upright, blocked/halved) by pattern
	x = np.arange(len(patterns))
	width = 0.2
	fig1, ax1 = plt.subplots(figsize=(12, 8))
	ax1.bar(x - 1.5*width, icp_supine_blocked, width, label='Supine Gx=0', color='black', edgecolor='black')
	ax1.bar(x - 0.5*width, icp_upright_blocked, width, label='Upright Gx=0', color='gray', edgecolor='black')
	ax1.bar(x + 0.5*width, icp_supine_halved, width, label='Supine Gx/2', color='white', edgecolor='black', hatch='//')
	ax1.bar(x + 1.5*width, icp_upright_halved, width, label='Upright Gx/2', color='white', edgecolor='black', hatch='..')
	ax1.set_xlabel('Pattern', fontsize=12, fontweight='bold')
	ax1.set_ylabel('ICP (mmHg)', fontsize=12, fontweight='bold')
	ax1.set_title('Stenotic Patterns Analysis: ICP', fontsize=14, fontweight='bold')
	ax1.set_xticks(x)
	ax1.set_xticklabels(patterns, fontsize=11)
	ax1.legend(fontsize=11, framealpha=0.9)
	ax1.grid(True, alpha=0.3, axis='y')
	plt.tight_layout()
	plt.savefig('Gadda_StenosisPatterns_ICP.png', dpi=300, bbox_inches='tight')
	print("\n" + "="*70)
	print("Plot saved as 'Gadda_StenosisPatterns_ICP.png'")
	print("="*70)
	plt.show(block=False)

	fig2, ax2 = plt.subplots(figsize=(12, 8))
	ax2.bar(x - 1.5*width, pvs_supine_blocked, width, label='Supine Gx=0', color='black', edgecolor='black')
	ax2.bar(x - 0.5*width, pvs_upright_blocked, width, label='Upright Gx=0', color='gray', edgecolor='black')
	ax2.bar(x + 0.5*width, pvs_supine_halved, width, label='Supine Gx/2', color='white', edgecolor='black', hatch='//')
	ax2.bar(x + 1.5*width, pvs_upright_halved, width, label='Upright Gx/2', color='white', edgecolor='black', hatch='..')
	ax2.set_xlabel('Pattern', fontsize=12, fontweight='bold')
	ax2.set_ylabel('Pvs (mmHg)', fontsize=12, fontweight='bold')
	ax2.set_title('Stenotic Patterns Analysis: Pvs', fontsize=14, fontweight='bold')
	ax2.set_xticks(x)
	ax2.set_xticklabels(patterns, fontsize=11)
	ax2.legend(fontsize=11, framealpha=0.9)
	ax2.grid(True, alpha=0.3, axis='y')
	plt.tight_layout()
	plt.savefig('Gadda_StenosisPatterns_Pvs.png', dpi=300, bbox_inches='tight')
	print("\n" + "="*70)
	print("Plot saved as 'Gadda_StenosisPatterns_Pvs.png'")
	print("="*70)
	plt.show(block=False)

	return (icp_supine_blocked, icp_upright_blocked, icp_supine_halved, icp_upright_halved,
			pvs_supine_blocked, pvs_upright_blocked, pvs_supine_halved, pvs_upright_halved)

# Comprehensive mass balance check
def comprehensive_mass_balance_check():
	"""
	Perform a comprehensive mass balance check for the model under various drainage pathway occlusion conditions.
	Prints inflow, outflow, and mass balance for each condition and posture (supine/upright).
	"""
	print("\n" + "="*70)
	print("COMPREHENSIVE MASS BALANCE CHECK")
	print("="*70)
	conditions = ['Basal', 'Gjr3=0', 'Gjr2=0', 'Gjr1=0', 'Gvvr=0', 'Gc3=0']
	for condition in conditions:
		print(f"\n{'='*70}")
		print(f"CONDITION: {condition}")
		print(f"{'='*70}")
		for posture in ['supine', 'upright']:
			print(f"\n{posture.upper()}:")
			print("-"*70)
			params = ModelParameters(posture=posture)
			if condition == 'Gjr3=0': params.kjr3 = 0.0; params.kjl3 = 0.0
			elif condition == 'Gjr2=0': params.kjr2 = 0.0; params.kjl2 = 0.0
			elif condition == 'Gjr1=0': params.kjr1 = 0.0; params.kjl1 = 0.0
			elif condition == 'Gvvr=0': params.Gvvr = 0.0; params.Gvvl = 0.0
			elif condition == 'Gc3=0': params.Gc3 = 0.0
			y0 = params.get_initial_conditions()
			sol = solve_ivp(gadda_ode_system_FIXED, [0, 500], y0, args=(params,),
						   method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8)
			if not sol.success:
				print("✗ Simulation failed!")
				continue
			final = sol.y[:, -1]
			Pic = final[0]; Ppa = final[1]; Pv = final[2]; Pvs = final[3]
			Pjr3 = final[4]; Pjl3 = final[5]; Pjr2 = final[6]; Pjl2 = final[7]
			Pc3 = final[8]; Pc2 = final[9]; Pvv = final[10]; Pazy = final[11]; Psvc = final[12]; Cpa = final[14]
			Rpa = cerebral_arterioles_resistance(Cpa, Ppa, Pic, params)
			Q_arterial = (params.Pa - Ppa) / (params.Rla + Rpa / 2.0)
			Pc_num = (Pv / params.Rpv) + (Ppa / (Rpa / 2.0)) + (Pic / params.Rf)
			Pc_den = (1.0 / params.Rpv) + (1.0 / (Rpa / 2.0)) + (1.0 / params.Rf)
			Pc = Pc_num / Pc_den
			Qf = (Pc - Pic) / params.Rf if Pc > Pic else 0.0
			Q0 = (Pic - Pvs) * params.G0 if Pic > Pvs else 0.0
			Qj3 = (Pvs - Pjr3) * jugular_conductance(Pjr3, params.Pj3ext, params.kjr3, params.A) + (Pvs - Pjl3) * jugular_conductance(Pjl3, params.Pj3ext, params.kjl3, params.A)
			Qvv_from_pvs = (Pvs - Pvv) * params.Gvvr + (Pvs - Pvv) * params.Gvvl
			Qc3_to_pvs = (Pvs - Pc3) * params.Gc3
			Qext = (params.Pa - Pc3) * params.Gex
			Q_out_direct = Qj3 + Qvv_from_pvs + Qc3_to_pvs
			Q_to_heart = Q_out_direct + Qext
			balance_pvs = Q_arterial - Q0 - Q_out_direct
			balance_heart = Q_arterial + Qext - Q_to_heart
			print(f"  INFLOWS:")
			print(f"    Q (cerebral arterial):        {Q_arterial:8.3f} ml/s")
			print(f"    Qext (external carotid):      {Qext:8.3f} ml/s")
			print(f"    Qf (CSF formation):           {Qf:8.3f} ml/s")
			print(f"    TOTAL IN:                     {Q_arterial + Qext:8.3f} ml/s")
			print(f"\n  OUTFLOWS FROM Pvs:")
			print(f"    Qj3 (to jugular J3):          {Qj3:8.3f} ml/s")
			print(f"    Qvv (to vertebral):           {Qvv_from_pvs:8.3f} ml/s")
			print(f"    Qc3 (to collateral):          {Qc3_to_pvs:8.3f} ml/s")
			print(f"    Q0 (CSF absorption):          {Q0:8.3f} ml/s")
			print(f"    TOTAL OUT from Pvs:           {Q_out_direct + Q0:8.3f} ml/s")
			print(f"\n  MASS BALANCE:")
			print(f"    Arterial IN - Venous OUT:     {balance_pvs:+8.3f} ml/s")
			print(f"    Total IN - Total to Heart:    {balance_heart:+8.3f} ml/s")
			print(f"    Pvs pressure:                 {Pvs:8.3f} mmHg")
			if abs(balance_heart) < 0.5:
				print(f"    ✓ MASS BALANCE OK")
			else:
				print(f"    ⚠ MASS BALANCE ISSUE!")
# Qvv (vertebral outflow) sensitivity analysis - BILATERAL
def sensitivity_analysis_qvv_bilateral():
	"""
	Sensitivity analysis of vertebral outflow (Qvv) under bilateral blockage conditions.
	Tests the effect of blocking each drainage pathway on both sides for Qvv in supine and upright postures.
	"""
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS: VERTEBRAL VEINS OUTFLOW (Qvv) - BILATERAL")
	print("="*70)
	conditions = ['Basal', 'Gj3=0', 'Gj2=0', 'Gj1=0', 'Gvv=0', 'Gc3=0']
	qvv_supine_values = []
	qvv_upright_values = []
	for condition in conditions:
		print(f"\n{'='*70}")
		print(f"Condition: {condition}")
		print(f"{'='*70}")
		# --- SUPINE ---
		print("\nSupine simulation...")
		params_supine = ModelParameters(posture='supine')
		if condition == 'Gj3=0':
			params_supine.kjr3 = 0.0
			params_supine.kjl3 = 0.0
		elif condition == 'Gj2=0':
			params_supine.kjr2 = 0.0
			params_supine.kjl2 = 0.0
		elif condition == 'Gj1=0':
			params_supine.kjr1 = 0.0
			params_supine.kjl1 = 0.0
		elif condition == 'Gvv=0':
			params_supine.Gvvr = 0.0
			params_supine.Gvvl = 0.0
		elif condition == 'Gc3=0':
			params_supine.Gc3 = 0.0
		y0 = params_supine.get_initial_conditions()
		sol_supine = solve_ivp(
			gadda_ode_system_FIXED, [0, 500], y0, args=(params_supine,),
			method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
		)
		if sol_supine.success:
			final = sol_supine.y[:, -1]
			Pvs_sup = final[3]
			Pvv_sup = final[10]
			Qvv_sup = (Pvs_sup - Pvv_sup) * params_supine.Gvvr + (Pvs_sup - Pvv_sup) * params_supine.Gvvl
			qvv_supine_values.append(Qvv_sup)
			print(f"✓ Supine Qvv = {Qvv_sup:.3f} ml/s")
		else:
			print("✗ Supine simulation failed!")
			qvv_supine_values.append(np.nan)
		# --- UPRIGHT ---
		print("Upright simulation...")
		params_upright = ModelParameters(posture='upright')
		if condition == 'Gj3=0':
			params_upright.kjr3 = 0.0
			params_upright.kjl3 = 0.0
		elif condition == 'Gj2=0':
			params_upright.kjr2 = 0.0
			params_upright.kjl2 = 0.0
		elif condition == 'Gj1=0':
			params_upright.kjr1 = 0.0
			params_upright.kjl1 = 0.0
		elif condition == 'Gvv=0':
			params_upright.Gvvr = 0.0
			params_upright.Gvvl = 0.0
		elif condition == 'Gc3=0':
			params_upright.Gc3 = 0.0
		y0 = params_upright.get_initial_conditions()
		sol_upright = solve_ivp(
			gadda_ode_system_FIXED, [0, 500], y0, args=(params_upright,),
			method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
		)
		if sol_upright.success:
			final = sol_upright.y[:, -1]
			Pvs_upr = final[3]
			Pvv_upr = final[10]
			Qvv_upr = (Pvs_upr - Pvv_upr) * params_upright.Gvvr + (Pvs_upr - Pvv_upr) * params_upright.Gvvl
			qvv_upright_values.append(Qvv_upr)
			print(f"✓ Upright Qvv = {Qvv_upr:.3f} ml/s")
		else:
			print("✗ Upright simulation failed!")
			qvv_upright_values.append(np.nan)
	# Print summary
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS SUMMARY")
	print("="*70)
	print(f"{'Condition':<15} {'Supine Qvv':<15} {'Upright Qvv':<15} {'Δ Qvv':<15}")
	print("-"*70)
	for i, cond in enumerate(conditions):
		delta = qvv_upright_values[i] - qvv_supine_values[i]
		print(f"{cond:<15} {qvv_supine_values[i]:>12.3f}   {qvv_upright_values[i]:>12.3f}   {delta:>+12.3f}")
	print("="*70)
	# Create bar plot
	fig, ax = plt.subplots(figsize=(12, 8))
	x = np.arange(len(conditions))
	width = 0.4
	bars1 = ax.bar(x - width/2, qvv_supine_values, width,
				   label='Supine', color='black', edgecolor='black', linewidth=1.5)
	bars2 = ax.bar(x + width/2, qvv_upright_values, width,
				   label='Upright', color='white', edgecolor='black', linewidth=1.5)
	ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
	ax.set_ylabel('Qvv (ml/s)', fontsize=12, fontweight='bold')
	ax.set_title('Sensitivity Analysis: Vertebral Veins Outflow (Qvv) - BILATERAL', 
				 fontsize=14, fontweight='bold')
	ax.set_xticks(x)
	ax.set_xticklabels(conditions, fontsize=11, rotation=0)
	ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
	ax.grid(True, alpha=0.3, axis='y')
	ax.set_ylim(0, None)
	plt.tight_layout()
	plt.savefig('Gadda_Qvv_sensitivity_bilateral.png', dpi=300, bbox_inches='tight')
	print("\n" + "="*70)
	print("Plot saved as 'Gadda_Qvv_sensitivity_bilateral.png'")
	print("="*70)
	plt.show(block=False)
	return qvv_supine_values, qvv_upright_values

# Qvv (vertebral outflow) sensitivity analysis - UNILATERAL
def sensitivity_analysis_qvv_unilateral():
	"""
	Sensitivity analysis of vertebral outflow (Qvv) under unilateral (right-side) blockage conditions.
	Tests the effect of right-side only occlusion for Qvv in supine and upright postures.
	"""
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS: VERTEBRAL VEINS OUTFLOW (Qvv) - UNILATERAL (Right Only)")
	print("="*70)
	conditions = ['Basal', 'Gjr3=0', 'Gjr2=0', 'Gjr1=0', 'Gvvr=0', 'Gc3=0']
	qvv_supine_values = []
	qvv_upright_values = []
	for condition in conditions:
		print(f"\n{'='*70}")
		print(f"Condition: {condition}")
		print(f"{'='*70}")
		# SUPINE
		print("\nSupine simulation...")
		params_sup = ModelParameters(posture='supine')
		if condition == 'Gjr3=0': params_sup.kjr3 = 0.0
		elif condition == 'Gjr2=0': params_sup.kjr2 = 0.0
		elif condition == 'Gjr1=0': params_sup.kjr1 = 0.0
		elif condition == 'Gvvr=0': params_sup.Gvvr = 0.0
		elif condition == 'Gc3=0': params_sup.Gc3 = 0.0
		y0 = params_sup.get_initial_conditions()
		sol_sup = solve_ivp(gadda_ode_system_FIXED, [0, 500], y0, args=(params_sup,),
						   method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8)
		if sol_sup.success:
			final = sol_sup.y[:, -1]
			Pvs, Pvv = final[3], final[10]
			Qvv = (Pvs - Pvv) * params_sup.Gvvr + (Pvs - Pvv) * params_sup.Gvvl
			qvv_supine_values.append(Qvv)
			print(f"✓ Supine Qvv = {Qvv:.3f} ml/s")
		else:
			qvv_supine_values.append(np.nan)
			print("✗ Supine failed!")
		# UPRIGHT
		print("Upright simulation...")
		params_upr = ModelParameters(posture='upright')
		if condition == 'Gjr3=0': params_upr.kjr3 = 0.0
		elif condition == 'Gjr2=0': params_upr.kjr2 = 0.0
		elif condition == 'Gjr1=0': params_upr.kjr1 = 0.0
		elif condition == 'Gvvr=0': params_upr.Gvvr = 0.0
		elif condition == 'Gc3=0': params_upr.Gc3 = 0.0
		y0 = params_upr.get_initial_conditions()
		sol_upr = solve_ivp(gadda_ode_system_FIXED, [0, 500], y0, args=(params_upr,),
						   method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8)
		if sol_upr.success:
			final = sol_upr.y[:, -1]
			Pvs, Pvv = final[3], final[10]
			Qvv = (Pvs - Pvv) * params_upr.Gvvr + (Pvs - Pvv) * params_upr.Gvvl
			qvv_upright_values.append(Qvv)
			print(f"✓ Upright Qvv = {Qvv:.3f} ml/s")
		else:
			qvv_upright_values.append(np.nan)
			print("✗ Upright failed!")
	# Summary
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS SUMMARY - UNILATERAL")
	print("="*70)
	print(f"{'Condition':<15} {'Supine Qvv':<15} {'Upright Qvv':<15} {'Δ Qvv':<15}")
	print("-"*70)
	for i, cond in enumerate(conditions):
		delta = qvv_upright_values[i] - qvv_supine_values[i]
		print(f"{cond:<15} {qvv_supine_values[i]:>12.3f}   {qvv_upright_values[i]:>12.3f}   {delta:>+12.3f}")
	print("="*70)
	# Plot
	fig, ax = plt.subplots(figsize=(12, 8))
	x = np.arange(len(conditions))
	width = 0.4
	ax.bar(x - width/2, qvv_supine_values, width, label='Supine', color='black', edgecolor='black', linewidth=1.5)
	ax.bar(x + width/2, qvv_upright_values, width, label='Upright', color='white', edgecolor='black', linewidth=1.5)
	ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
	ax.set_ylabel('Qvv (ml/s)', fontsize=12, fontweight='bold')
	ax.set_title('Sensitivity Analysis: Vertebral Veins Outflow (Qvv) - UNILATERAL (Right)', fontsize=14, fontweight='bold')
	ax.set_xticks(x)
	ax.set_xticklabels(conditions, fontsize=11, rotation=0)
	ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
	ax.grid(True, alpha=0.3, axis='y')
	ax.set_ylim(0, None)
	plt.tight_layout()
	plt.savefig('Gadda_Qvv_sensitivity_unilateral.png', dpi=300, bbox_inches='tight')
	print("\n" + "="*70)
	print("Plot saved as 'Gadda_Qvv_sensitivity_unilateral.png'")
	print("="*70)
	plt.show(block=False)
	return qvv_supine_values, qvv_upright_values
# Qsvc1 (jugular outflow) sensitivity analysis - UNILATERAL
def sensitivity_analysis_qsvc_unilateral():
	"""
	Sensitivity analysis of jugular outflow (Qsvc) under unilateral (right-side) blockage conditions.
	Tests the effect of right-side only occlusion for Qsvc in supine and upright postures.
	"""
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS: JUGULAR VEINS OUTFLOW (Qsvc1) - UNILATERAL (Right Only)")
	print("="*70)
	conditions = ['Basal', 'Gjr3=0', 'Gjr2=0', 'Gjr1=0', 'Gvvr=0', 'Gc3=0']
	qsvc_supine_values = []
	qsvc_upright_values = []
	for condition in conditions:
		print(f"\n{'='*70}")
		print(f"Condition: {condition}")
		print(f"{'='*70}")
		# SUPINE
		print("\nSupine simulation...")
		params_sup = ModelParameters(posture='supine')
		if condition == 'Gjr3=0': params_sup.kjr3 = 0.0
		elif condition == 'Gjr2=0': params_sup.kjr2 = 0.0
		elif condition == 'Gjr1=0': params_sup.kjr1 = 0.0
		elif condition == 'Gvvr=0': params_sup.Gvvr = 0.0
		elif condition == 'Gc3=0': params_sup.Gc3 = 0.0
		y0 = params_sup.get_initial_conditions()
		sol_sup = solve_ivp(gadda_ode_system_FIXED, [0, 500], y0, args=(params_sup,),
						   method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8)
		if sol_sup.success:
			final = sol_sup.y[:, -1]
			Pjr2, Pjl2, Psvc = final[6], final[7], final[12]
			Gjr1 = jugular_conductance(Pjr2, params_sup.Pj1ext, params_sup.kjr1, params_sup.A)
			Gjl1 = jugular_conductance(Pjl2, params_sup.Pj1ext, params_sup.kjl1, params_sup.A)
			Qsvc1 = (Pjr2 - Psvc) * Gjr1 + (Pjl2 - Psvc) * Gjl1
			qsvc_supine_values.append(Qsvc1)
			print(f"✓ Supine Qsvc1 = {Qsvc1:.3f} ml/s")
		else:
			qsvc_supine_values.append(np.nan)
			print("✗ Supine failed!")
		# UPRIGHT
		print("Upright simulation...")
		params_upr = ModelParameters(posture='upright')
		if condition == 'Gjr3=0': params_upr.kjr3 = 0.0
		elif condition == 'Gjr2=0': params_upr.kjr2 = 0.0
		elif condition == 'Gjr1=0': params_upr.kjr1 = 0.0
		elif condition == 'Gvvr=0': params_upr.Gvvr = 0.0
		elif condition == 'Gc3=0': params_upr.Gc3 = 0.0
		y0 = params_upr.get_initial_conditions()
		sol_upr = solve_ivp(gadda_ode_system_FIXED, [0, 500], y0, args=(params_upr,),
						   method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8)
		if sol_upr.success:
			final = sol_upr.y[:, -1]
			Pjr2, Pjl2, Psvc = final[6], final[7], final[12]
			Gjr1 = jugular_conductance(Pjr2, params_upr.Pj1ext, params_upr.kjr1, params_upr.A)
			Gjl1 = jugular_conductance(Pjl2, params_upr.Pj1ext, params_upr.kjl1, params_upr.A)
			Qsvc1 = (Pjr2 - Psvc) * Gjr1 + (Pjl2 - Psvc) * Gjl1
			qsvc_upright_values.append(Qsvc1)
			print(f"✓ Upright Qsvc1 = {Qsvc1:.3f} ml/s")
		else:
			qsvc_upright_values.append(np.nan)
			print("✗ Upright failed!")
	# Summary
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS SUMMARY - UNILATERAL")
	print("="*70)
	print(f"{'Condition':<15} {'Supine Qsvc1':<15} {'Upright Qsvc1':<15} {'Δ Qsvc1':<15}")
	print("-"*70)
	for i, cond in enumerate(conditions):
		delta = qsvc_upright_values[i] - qsvc_supine_values[i]
		print(f"{cond:<15} {qsvc_supine_values[i]:>12.3f}   {qsvc_upright_values[i]:>12.3f}   {delta:>+12.3f}")
	print("="*70)
	# Plot
	fig, ax = plt.subplots(figsize=(12, 8))
	x = np.arange(len(conditions))
	width = 0.4
	ax.bar(x - width/2, qsvc_supine_values, width, label='Supine', color='black', edgecolor='black', linewidth=1.5)
	ax.bar(x + width/2, qsvc_upright_values, width, label='Upright', color='white', edgecolor='black', linewidth=1.5)
	ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
	ax.set_ylabel('Qsvc1 (ml/s)', fontsize=12, fontweight='bold')
	ax.set_title('Sensitivity Analysis: Jugular Veins Outflow (Qsvc1) - UNILATERAL (Right)', fontsize=14, fontweight='bold')
	ax.set_xticks(x)
	ax.set_xticklabels(conditions, fontsize=11, rotation=0)
	ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
	ax.grid(True, alpha=0.3, axis='y')
	ax.set_ylim(0, None)
	plt.tight_layout()
	plt.savefig('Gadda_Qsvc1_sensitivity_unilateral.png', dpi=300, bbox_inches='tight')
	print("\n" + "="*70)
	print("Plot saved as 'Gadda_Qsvc1_sensitivity_unilateral.png'")
	print("="*70)
	plt.show(block=False)
	return qsvc_supine_values, qsvc_upright_values
# Qsvc1 (jugular outflow) sensitivity analysis - BILATERAL
def sensitivity_analysis_qsvc_bilateral():
	"""
	Sensitivity analysis of jugular outflow (Qsvc) under bilateral blockage conditions.
	Tests the effect of blocking each drainage pathway on both sides for Qsvc in supine and upright postures.
	"""
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS: JUGULAR VEINS OUTFLOW (Qsvc1) - BILATERAL")
	print("="*70)
	conditions = ['Basal', 'Gj3=0', 'Gj2=0', 'Gj1=0', 'Gvv=0', 'Gc3=0']
	qsvc_supine_values = []
	qsvc_upright_values = []
	for condition in conditions:
		print(f"\n{'='*70}")
		print(f"Condition: {condition}")
		print(f"{'='*70}")
		# --- SUPINE ---
		print("\nSupine simulation...")
		params_supine = ModelParameters(posture='supine')
		if condition == 'Gj3=0':
			params_supine.kjr3 = 0.0
			params_supine.kjl3 = 0.0
		elif condition == 'Gj2=0':
			params_supine.kjr2 = 0.0
			params_supine.kjl2 = 0.0
		elif condition == 'Gj1=0':
			params_supine.kjr1 = 0.0
			params_supine.kjl1 = 0.0
		elif condition == 'Gvv=0':
			params_supine.Gvvr = 0.0
			params_supine.Gvvl = 0.0
		elif condition == 'Gc3=0':
			params_supine.Gc3 = 0.0
		y0 = params_supine.get_initial_conditions()
		sol_supine = solve_ivp(
			gadda_ode_system_FIXED, [0, 500], y0, args=(params_supine,),
			method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
		)
		if sol_supine.success:
			final = sol_supine.y[:, -1]
			Pjr2_sup = final[6]
			Pjl2_sup = final[7]
			Psvc_sup = final[12]
			# Jugular conductances
			Gjr1_sup = jugular_conductance(Pjr2_sup, params_supine.Pj1ext, params_supine.kjr1, params_supine.A)
			Gjl1_sup = jugular_conductance(Pjl2_sup, params_supine.Pj1ext, params_supine.kjl1, params_supine.A)
			Qsvc1_sup = (Pjr2_sup - Psvc_sup) * Gjr1_sup + (Pjl2_sup - Psvc_sup) * Gjl1_sup
			qsvc_supine_values.append(Qsvc1_sup)
			print(f"✓ Supine Qsvc1 = {Qsvc1_sup:.3f} ml/s")
		else:
			print("✗ Supine simulation failed!")
			qsvc_supine_values.append(np.nan)
		# --- UPRIGHT ---
		print("Upright simulation...")
		params_upright = ModelParameters(posture='upright')
		if condition == 'Gj3=0':
			params_upright.kjr3 = 0.0
			params_upright.kjl3 = 0.0
		elif condition == 'Gj2=0':
			params_upright.kjr2 = 0.0
			params_upright.kjl2 = 0.0
		elif condition == 'Gj1=0':
			params_upright.kjr1 = 0.0
			params_upright.kjl1 = 0.0
		elif condition == 'Gvv=0':
			params_upright.Gvvr = 0.0
			params_upright.Gvvl = 0.0
		elif condition == 'Gc3=0':
			params_upright.Gc3 = 0.0
		y0 = params_upright.get_initial_conditions()
		sol_upright = solve_ivp(
			gadda_ode_system_FIXED, [0, 500], y0, args=(params_upright,),
			method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
		)
		if sol_upright.success:
			final = sol_upright.y[:, -1]
			Pjr2_upr = final[6]
			Pjl2_upr = final[7]
			Psvc_upr = final[12]
			Gjr1_upr = jugular_conductance(Pjr2_upr, params_upright.Pj1ext, params_upright.kjr1, params_upright.A)
			Gjl1_upr = jugular_conductance(Pjl2_upr, params_upright.Pj1ext, params_upright.kjl1, params_upright.A)
			Qsvc1_upr = (Pjr2_upr - Psvc_upr) * Gjr1_upr + (Pjl2_upr - Psvc_upr) * Gjl1_upr
			qsvc_upright_values.append(Qsvc1_upr)
			print(f"✓ Upright Qsvc1 = {Qsvc1_upr:.3f} ml/s")
		else:
			print("✗ Upright simulation failed!")
			qsvc_upright_values.append(np.nan)
	# Print summary
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS SUMMARY")
	print("="*70)
	print(f"{'Condition':<15} {'Supine Qsvc1':<15} {'Upright Qsvc1':<15} {'Δ Qsvc1':<15}")
	print("-"*70)
	for i, cond in enumerate(conditions):
		delta = qsvc_upright_values[i] - qsvc_supine_values[i]
		print(f"{cond:<15} {qsvc_supine_values[i]:>12.3f}   {qsvc_upright_values[i]:>12.3f}   {delta:>+12.3f}")
	print("="*70)
	# Create bar plot
	fig, ax = plt.subplots(figsize=(12, 8))
	x = np.arange(len(conditions))
	width = 0.4
	bars1 = ax.bar(x - width/2, qsvc_supine_values, width,
				   label='Supine', color='black', edgecolor='black', linewidth=1.5)
	bars2 = ax.bar(x + width/2, qsvc_upright_values, width,
				   label='Upright', color='white', edgecolor='black', linewidth=1.5)
	ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
	ax.set_ylabel('Qsvc1 (ml/s)', fontsize=12, fontweight='bold')
	ax.set_title('Sensitivity Analysis: Jugular Veins Outflow (Qsvc1) - BILATERAL', 
				 fontsize=14, fontweight='bold')
	ax.set_xticks(x)
	ax.set_xticklabels(conditions, fontsize=11, rotation=0)
	ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
	ax.grid(True, alpha=0.3, axis='y')
	ax.set_ylim(0, None)
	plt.tight_layout()
	plt.savefig('Gadda_Qsvc1_sensitivity_bilateral.png', dpi=300, bbox_inches='tight')
	print("\n" + "="*70)
	print("Plot saved as 'Gadda_Qsvc1_sensitivity_bilateral.png'")
	print("="*70)
	plt.show(block=False)
	return qsvc_supine_values, qsvc_upright_values
def sensitivity_analysis_pvs_unilateral():
	"""
	Sensitivity analysis of venous sinus pressure (Pvs) under unilateral (right-side) blockage conditions.
	Tests the effect of right-side only occlusion for Pvs in supine and upright postures.
	"""
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS: VENOUS SINUSES PRESSURE (Pvs) - UNILATERAL (Right Only)")
	print("="*70)
	conditions = ['Basal', 'Gjr3=0', 'Gjr2=0', 'Gjr1=0', 'Gvvr=0', 'Gc3=0']
	pvs_supine_values = []
	pvs_upright_values = []
	for condition in conditions:
		print(f"\n{'='*70}")
		print(f"Condition: {condition}")
		print(f"{'='*70}")
		# --- SUPINE ---
		print("\nSupine simulation...")
		params_supine = ModelParameters(posture='supine')
		# Apply condition (UNILATERAL - RIGHT side only)
		if condition == 'Gjr3=0':
			params_supine.kjr3 = 0.0
		elif condition == 'Gjr2=0':
			params_supine.kjr2 = 0.0
		elif condition == 'Gjr1=0':
			params_supine.kjr1 = 0.0
		elif condition == 'Gvvr=0':
			params_supine.Gvvr = 0.0
		elif condition == 'Gc3=0':
			params_supine.Gc3 = 0.0
		y0 = params_supine.get_initial_conditions()
		sol_supine = solve_ivp(
			gadda_ode_system_FIXED, [0, 500], y0, args=(params_supine,),
			method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
		)
		if sol_supine.success:
			pvs_sup = sol_supine.y[3, -1]
			pvs_supine_values.append(pvs_sup)
			print(f"✓ Supine Pvs = {pvs_sup:.3f} mmHg")
		else:
			print("✗ Supine simulation failed!")
			pvs_supine_values.append(np.nan)
		# --- UPRIGHT ---
		print("Upright simulation...")
		params_upright = ModelParameters(posture='upright')
		if condition == 'Gjr3=0':
			params_upright.kjr3 = 0.0
		elif condition == 'Gjr2=0':
			params_upright.kjr2 = 0.0
		elif condition == 'Gjr1=0':
			params_upright.kjr1 = 0.0
		elif condition == 'Gvvr=0':
			params_upright.Gvvr = 0.0
		elif condition == 'Gc3=0':
			params_upright.Gc3 = 0.0
		y0 = params_upright.get_initial_conditions()
		sol_upright = solve_ivp(
			gadda_ode_system_FIXED, [0, 500], y0, args=(params_upright,),
			method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
		)
		if sol_upright.success:
			pvs_upr = sol_upright.y[3, -1]
			pvs_upright_values.append(pvs_upr)
			print(f"✓ Upright Pvs = {pvs_upr:.3f} mmHg")
		else:
			print("✗ Upright simulation failed!")
			pvs_upright_values.append(np.nan)
	# Print summary
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS SUMMARY - UNILATERAL")
	print("="*70)
	print(f"{'Condition':<15} {'Supine Pvs':<15} {'Upright Pvs':<15} {'Δ Pvs':<15}")
	print("-"*70)
	for i, cond in enumerate(conditions):
		delta = pvs_upright_values[i] - pvs_supine_values[i]
		print(f"{cond:<15} {pvs_supine_values[i]:>12.3f}   {pvs_upright_values[i]:>12.3f}   {delta:>+12.3f}")
	print("="*70)
	# Create bar plot
	fig, ax = plt.subplots(figsize=(12, 8))
	x = np.arange(len(conditions))
	width = 0.4
	bars1 = ax.bar(x - width/2, pvs_supine_values, width, label='Supine', color='black', edgecolor='black', linewidth=1.5)
	bars2 = ax.bar(x + width/2, pvs_upright_values, width, label='Upright', color='white', edgecolor='black', linewidth=1.5)
	ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
	ax.set_ylabel('Pvs (mmHg)', fontsize=12, fontweight='bold')
	ax.set_title('Sensitivity Analysis: Venous Sinuses Pressure (Pvs) - UNILATERAL (Right)', fontsize=14, fontweight='bold')
	ax.set_xticks(x)
	ax.set_xticklabels(conditions, fontsize=11, rotation=0)
	ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
	ax.grid(True, alpha=0.3, axis='y')
	ax.set_ylim(0, None)
	plt.tight_layout()
	plt.savefig('Gadda_Pvs_sensitivity_unilateral.png', dpi=300, bbox_inches='tight')
	print("\n" + "="*70)
	print("Plot saved as 'Gadda_Pvs_sensitivity_unilateral.png'")
	print("="*70)
	plt.show(block=False)
	return pvs_supine_values, pvs_upright_values

"""
Full Gadda et al. (2015) model reproduction using modular codebase.
Reproduces all baseline, posture, sensitivity, stenosis, and mass balance analyses for publication.
"""

import numpy as np
import matplotlib.pyplot as plt
from gadda_model.equations import (
	gadda_ode_system_FIXED,
	jugular_conductance,
	cerebral_arterioles_resistance,
)
from gadda_model.parameters import ModelParameters
from scipy.integrate import solve_ivp

def run_baseline(posture='supine', t_span=[0, 300]):
	"""
	Run a baseline simulation for a given posture ('supine' or 'upright') over a specified time span.
	Returns the ODE solution and model parameters.
	"""
	params = ModelParameters(posture=posture)
	y0 = params.get_initial_conditions()
	sol = solve_ivp(
		gadda_ode_system_FIXED, t_span, y0, args=(params,),
		method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
	)
	return sol, params

def print_baseline_results(sol, params, posture):
	"""
	Print summary results for a baseline simulation, including ICP, Pvs, jugular and vertebral outflows, and total cerebral blood flow.
	"""
	final = sol.y[:, -1]
	Pic, Pvs, Pjr3, Pjl3, Pvv = final[0], final[3], final[4], final[5], final[10]
	Cpa = final[14]
	Ppa = final[1]
	Rpa = cerebral_arterioles_resistance(Cpa, Ppa, Pic, params)
	Q_cerebral = (params.Pa - Ppa) / (params.Rla + Rpa / 2.0)
	# Jugular conductances
	Gjr3 = jugular_conductance(Pjr3, params.Pj3ext, params.kjr3, params.A)
	Gjl3 = jugular_conductance(Pjl3, params.Pj3ext, params.kjl3, params.A)
	Qj3 = (Pvs - Pjr3) * Gjr3 + (Pvs - Pjl3) * Gjl3
	Qvv = (Pvs - Pvv) * params.Gvvr + (Pvs - Pvv) * params.Gvvl
	print(f"\nBaseline ({posture}):")
	print(f"  ICP: {Pic:.2f} mmHg")
	print(f"  Pvs: {Pvs:.2f} mmHg")
	print(f"  Qj3: {Qj3:.2f} ml/s")
	print(f"  Qvv: {Qvv:.2f} ml/s")
	print(f"  Q_cerebral: {Q_cerebral:.2f} ml/s")

def compare_supine_upright():
	"""
	Compare baseline model results for supine and upright postures.
	Prints summary statistics and generates a comparison plot of key flow variables.
	"""
	sol_sup, params_sup = run_baseline('supine')
	sol_upr, params_upr = run_baseline('upright')
	print("\n--- SUPINE ---")
	print_baseline_results(sol_sup, params_sup, 'supine')
	print("\n--- UPRIGHT ---")
	print_baseline_results(sol_upr, params_upr, 'upright')
	# Plot comparison (like Gadda Fig 4)
	# ...existing code for plotting comparison...

def simulate_posture_change(t_change=300, t_total=600):
	"""
	Simulate a dynamic posture change from supine to upright at a specified time.
	Plots and saves the time-course of ICP and Pvs across the transition.
	Returns the ODE solutions for both phases.
	"""
	params_sup = ModelParameters('supine')
	y0 = params_sup.get_initial_conditions()
	sol_sup = solve_ivp(
		gadda_ode_system_FIXED, [0, t_change], y0, args=(params_sup,),
		method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
	)
	params_upr = ModelParameters('upright')
	y0_upr = sol_sup.y[:, -1]
	sol_upr = solve_ivp(
		gadda_ode_system_FIXED, [t_change, t_total], y0_upr, args=(params_upr,),
		method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
	)

	# Combine results for plotting
	t_combined = np.concatenate([sol_sup.t, sol_upr.t])
	ICP_combined = np.concatenate([sol_sup.y[0], sol_upr.y[0]])
	Pvs_combined = np.concatenate([sol_sup.y[3], sol_upr.y[3]])

	# Plot ICP and Pvs over time with posture change
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

	ax1.plot(t_combined, ICP_combined, 'b-', linewidth=2)
	ax1.axvline(t_change, color='red', linestyle='--', linewidth=2, label='Posture Change')
	ax1.set_ylabel('ICP (mmHg)', fontsize=12, fontweight='bold')
	ax1.set_title('Posture Change: Supine → Upright', fontsize=14, fontweight='bold')
	ax1.grid(True, alpha=0.3)
	ax1.legend()

	ax2.plot(t_combined, Pvs_combined, 'g-', linewidth=2)
	ax2.axvline(t_change, color='red', linestyle='--', linewidth=2, label='Posture Change')
	ax2.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
	ax2.set_ylabel('Pvs (mmHg)', fontsize=12, fontweight='bold')
	ax2.grid(True, alpha=0.3)
	ax2.legend()

	plt.tight_layout()
	# Save and show plot
	import os
	os.makedirs('validation_plots', exist_ok=True)
	filename = f'validation_plots/validation_posture_change_{t_change}s.png'
	plt.savefig(filename, dpi=300, bbox_inches='tight')
	print(f"\nPlot saved: {filename}")
	plt.show(block=False)

	return sol_sup, sol_upr


def sensitivity_analysis_pvs_bilateral():
	"""
	Sensitivity analysis of venous sinus pressure (Pvs) under bilateral blockage conditions.
	Tests the effect of blocking each drainage pathway on both sides for Pvs in supine and upright postures.
	"""
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS: VENOUS SINUSES PRESSURE (Pvs) - BILATERAL")
	print("="*70)
	conditions = ['Basal', 'Gj3=0', 'Gj2=0', 'Gj1=0', 'Gvv=0', 'Gc3=0']
	pvs_supine_values = []
	pvs_upright_values = []
	for condition in conditions:
		print(f"\n{'='*70}")
		print(f"Condition: {condition}")
		print(f"{'='*70}")
		# --- SUPINE ---
		print("\nSupine simulation...")
		params_supine = ModelParameters(posture='supine')
		# Apply condition (BILATERAL - both sides)
		if condition == 'Gj3=0':
			params_supine.kjr3 = 0.0
			params_supine.kjl3 = 0.0
		elif condition == 'Gj2=0':
			params_supine.kjr2 = 0.0
			params_supine.kjl2 = 0.0
		elif condition == 'Gj1=0':
			params_supine.kjr1 = 0.0
			params_supine.kjl1 = 0.0
		elif condition == 'Gvv=0':
			params_supine.Gvvr = 0.0
			params_supine.Gvvl = 0.0
		elif condition == 'Gc3=0':
			params_supine.Gc3 = 0.0
		y0 = params_supine.get_initial_conditions()
		sol_supine = solve_ivp(
			gadda_ode_system_FIXED, [0, 500], y0, args=(params_supine,),
			method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
		)
		if sol_supine.success:
			pvs_sup = sol_supine.y[3, -1]
			pvs_supine_values.append(pvs_sup)
			print(f"✓ Supine Pvs = {pvs_sup:.3f} mmHg")
		else:
			print("✗ Supine simulation failed!")
			pvs_supine_values.append(np.nan)
		# --- UPRIGHT ---
		print("Upright simulation...")
		params_upright = ModelParameters(posture='upright')
		if condition == 'Gj3=0':
			params_upright.kjr3 = 0.0
			params_upright.kjl3 = 0.0
		elif condition == 'Gj2=0':
			params_upright.kjr2 = 0.0
			params_upright.kjl2 = 0.0
		elif condition == 'Gj1=0':
			params_upright.kjr1 = 0.0
			params_upright.kjl1 = 0.0
		elif condition == 'Gvv=0':
			params_upright.Gvvr = 0.0
			params_upright.Gvvl = 0.0
		elif condition == 'Gc3=0':
			params_upright.Gc3 = 0.0
		y0 = params_upright.get_initial_conditions()
		sol_upright = solve_ivp(
			gadda_ode_system_FIXED, [0, 500], y0, args=(params_upright,),
			method='RK45', dense_output=True, max_step=1.0, rtol=1e-6, atol=1e-8
		)
		if sol_upright.success:
			pvs_upr = sol_upright.y[3, -1]
			pvs_upright_values.append(pvs_upr)
			print(f"✓ Upright Pvs = {pvs_upr:.3f} mmHg")
		else:
			print("✗ Upright simulation failed!")
			pvs_upright_values.append(np.nan)
	# Print summary
	print("\n" + "="*70)
	print("SENSITIVITY ANALYSIS SUMMARY")
	print("="*70)
	print(f"{'Condition':<15} {'Supine Pvs':<15} {'Upright Pvs':<15} {'Δ Pvs':<15}")
	print("-"*70)
	for i, cond in enumerate(conditions):
		delta = pvs_upright_values[i] - pvs_supine_values[i]
		print(f"{cond:<15} {pvs_supine_values[i]:>12.3f}   {pvs_upright_values[i]:>12.3f}   {delta:>+12.3f}")
	print("="*70)
	# Create bar plot
	fig, ax = plt.subplots(figsize=(12, 8))
	x = np.arange(len(conditions))
	width = 0.4
	bars1 = ax.bar(x - width/2, pvs_supine_values, width,
				   label='Supine', color='black', edgecolor='black', linewidth=1.5)
	bars2 = ax.bar(x + width/2, pvs_upright_values, width,
				   label='Upright', color='white', edgecolor='black', linewidth=1.5)
	ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
	ax.set_ylabel('Pvs (mmHg)', fontsize=12, fontweight='bold')
	ax.set_title('Sensitivity Analysis: Venous Sinuses Pressure (Pvs) - BILATERAL', 
				 fontsize=14, fontweight='bold')
	ax.set_xticks(x)
	ax.set_xticklabels(conditions, fontsize=11, rotation=0)
	ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
	ax.grid(True, alpha=0.3, axis='y')
	ax.set_ylim(0, None)
	plt.tight_layout()
	plt.savefig('Gadda_Pvs_sensitivity_bilateral.png', dpi=300, bbox_inches='tight')
	print("\n" + "="*70)
	print("Plot saved as 'Gadda_Pvs_sensitivity_bilateral.png'")
	print("="*70)
	plt.show(block=False)
	return pvs_supine_values, pvs_upright_values

if __name__ == "__main__":
	print("\n==============================")
	print("GADDA MODEL FULL REPRODUCTION")
	print("==============================")
	compare_supine_upright()
	print("\n[Posture Change Simulation]")
	simulate_posture_change(t_change=2, t_total=10)
	print("\n[Sensitivity Analyses]")
	sensitivity_analysis_pvs_bilateral()
	sensitivity_analysis_pvs_unilateral()
	sensitivity_analysis_qsvc_bilateral()
	sensitivity_analysis_qsvc_unilateral()
	sensitivity_analysis_qvv_bilateral()
	sensitivity_analysis_qvv_unilateral()

	print("\n[Stenosis Patterns Analysis]")
	stenosis_patterns_pvs_analysis()

	print("\n[Comprehensive Mass Balance Check]")
	comprehensive_mass_balance_check()
