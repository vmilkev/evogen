import sys

sys.path.append('/media/vimi/958ab6f9-af1c-4df2-ba22-4788a76e7c79/vimi/PROJECTS/codebase/evo/release')

import evogen

def test_model_3():

        # one trait setup
        # EXAMPLE: 11.2. p.183.
        # y = Xb + Za + e

        # corr: [a], identity matrix
        
        model = evogen.Model()
        solver = evogen.Pcg()
        
        iR = [0.0041]

        iG1 = [0.1004]

        model.append_residual(iR, 1)

        model.append_observation("y.dat") #obs := 0

        model.append_effect("x.dat", 0)   #eff := 0
        model.append_effect("z.dat", 0.0) #eff := 1

        corr_eff = [1]

        eff_trate = [0, 1]
        obs_trate = 0

        model.append_corrstruct(iG1, 1, "corZ.dat", corr_eff)

        model.append_traitstruct(obs_trate, eff_trate)

        solver.append_model(model)

        #model.print()

        #----------------------------------

        solver.solve()
        
        solver.get_solution( "solution_model_3.dat" )

        model.clear()


if __name__ == '__main__':      
#def main():
        test_model_3()

