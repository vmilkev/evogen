import sys

sys.path.append('/media/vimi/958ab6f9-af1c-4df2-ba22-4788a76e7c79/vimi/PROJECTS/codebase/evo/release')

import evogen

def test_model_2():

        # one trait setup
        # EXAMPLE: 7.1. p.110.
        # y = Xb + Ua + Wm + Sp + e
        # corr1: [a m]
        # corr2: [p]

        # changes to model interpretation:
        # X moves to random effects but without covariance structure
        # corr1: [a m]
        # corr2: [p]
        # corr0 == 0: [b], because it is a fixed effect

        model = evogen.Model()
        solver = evogen.Pcg()

        iR = [0.0029]
        
        iG2 = [0.0076, 0.0034, 0.0034, 0.0126]

        iG3 = [0.025]
        
        model.append_residual(iR, 1);
        
        model.append_observation("y.dat") #obs := 0

        model.append_effect("x.dat", 0)   #eff := 0
        model.append_effect("s.dat", 0)   #eff := 1
        model.append_effect("w.dat", 0)   #eff := 2
        model.append_effect("u.dat", 0)   #eff := 3

        corr_eff_1 = [3, 2]
        corr_eff_2 = [1]

        eff_trate = [0, 3, 2, 1]
        obs_trate = 0;

        model.append_corrstruct(iG2, 2, "iA.dat", corr_eff_1)
        model.append_corrstruct(iG3, 1, "corS.dat", corr_eff_2)

        model.append_traitstruct(obs_trate, eff_trate)

        solver.append_model(model)

        #----------------------------------

        solver.solve()

        solver.get_solution( "solution_model_2.dat" )

        model.clear()


if __name__ == '__main__':      
#def main():
        test_model_2()


