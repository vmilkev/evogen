import sys

sys.path.append('release')

import evogen

def test_model_4():

	model = evogen.Model()
	solver = evogen.Pcg()

	iR = [0.041]
	iG1 = [0.001004]
        
	model.append_residual(iR, 1);

	model.append_observation("tests/data/model_4/obs_10000.dat"); # obs := 0

	model.append_effect("tests/data/model_4/obs_10000_snp_1000.txt", 0); # eff := 0

	model.append_effect("tests/data/model_4/fixed_10000.dat", 0.0);           # eff := 1
        
	corr_eff = [0]

	eff_trate = [1, 0]
	obs_trate = 0
        
	model.append_corrstruct(iG1, 1, "I", 1000, corr_eff)

	model.append_traitstruct(obs_trate, eff_trate)

	solver.append_model(model)

	solver.solve(3)
        
	solver.get_solution( "py_out_solution_model_4.dat" )

	model.clear()
	
if __name__ == '__main__':	
#def main():
	test_model_4()

