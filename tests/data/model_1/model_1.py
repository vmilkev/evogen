import sys

sys.path.append('/media/vimi/958ab6f9-af1c-4df2-ba22-4788a76e7c79/vimi/PROJECTS/codebase/evo/release')

import evogen

def test_model_1():

        # two-trait model
        # EXAMPLE: 5.1. p.72.
        # y1 = Xb1 + Za1 + e1
        # y1 = Xb2 + Za2 + e2
        
        model = evogen.Model()
        solver = evogen.Pcg()

        dummy_type1 = int(1)
        dummy_type2 = float(1.0)
        
        model.append_effect("z1.dat", dummy_type1)   #eff_0
        model.append_effect("z2.dat", dummy_type1)   #eff_1
        model.append_effect("x1.dat", dummy_type2)   #eff_2
        model.append_effect("x2.dat", dummy_type2)   #eff_3
        
        model.append_observation("y1.dat")           #obs_0
        model.append_observation("y2.dat")           #obs_1

        model.append_residual("iR.dat")

        corr_eff = [0, 1]

        model.append_corrstruct("iG.dat", "ainv.dat", corr_eff);

        eff_trate_1 = [ 2, 0 ]                       #trait 1
        obs_trate_1 = 0                              #trait 1

        eff_trate_2 = [ 3, 1 ]                       #trait 2
        obs_trate_2 = 1                              #trait 2

        model.append_traitstruct( obs_trate_1, eff_trate_1 )
        model.append_traitstruct( obs_trate_2, eff_trate_2 )

        solver.append_model( model )

        solver.solve()

        solver.get_solution( "solution_model_1.dat" )

        model.clear()


if __name__ == '__main__':      
#def main():
        test_model_1()

