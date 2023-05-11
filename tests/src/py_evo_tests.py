import sys

sys.path.append('release')

import evogen
import unittest

class TestEVO( unittest.TestCase ):

    def testModel( self ):

        model = evogen.Model()

        self.assertEqual( model.size_of("res"), 0 )
        self.assertEqual( model.size_of("obs"), 0 )
        self.assertEqual( model.size_of("eff"), 0 )

        iR = [0.0278, -0.0102, -0.0102, 0.0371]

        y1 = [ 4.5, 2.9, 3.9, 3.5, 5.0 ]
        y2 = [ 6.8, 5.0, 6.8, 6.0, 7.5 ]

        x1 = [ 1, 0, 0, 1, 0, 1, 1, 0, 1, 0 ]
        x2 = [ 1, 0, 0, 1, 0, 1, 1, 0, 1, 0 ]
        z1 = [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 ]
        z2 = [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 ]

        iA = [ 1.833, 0.5, 0, -0.667, 0, -1, 0, 0,
                0.5, 2, 0.5, 0, -1, -1, 0, 0,
                0, 0.5, 2, 0, -1, 0.5, 0, -1,
                -0.667, 0, 0, 1.833, 0.5, 0, -1, 0,
                0, -1, -1, 0.5, 2.5, 0, -1, 0,
                -1, -1, 0.5, 0, 0, 2.5, 0, -1,
                0, 0, 0, -1, -1, 0, 2, 0,
                0, 0, -1, 0, 0, -1, 0, 2 ]

        iG1 = [ 0.084, -0.0378, -0.0378, 0.042 ]
        corr_eff = [ 2, 3 ]

        eff_trate_1 = [ 0, 2 ]
        obs_trate_1 = 0

        eff_trate_2 = [ 1, 3 ]
        obs_trate_2 = 1

        model.append_residual( iR, 2)
        model.append_observation( y1, 5 )

        self.assertEqual( model.size_of("res"), 4 )
        self.assertEqual( model.size_of("obs"), 5 )

        model.append_residual( iR, 2)
        model.append_observation( y2, 5 )

        self.assertEqual( model.size_of("res"), 8 )
        self.assertEqual( model.size_of("obs"), 10 )

        model.append_effect( x1, 5, 2 )              #eff := 0
        self.assertEqual( model.size_of("eff"), 10 )

        model.append_effect( x2, 5, 2 )              #eff := 1
        self.assertEqual( model.size_of("eff"), 20 )

        model.append_effect( z1, 5, 8 )              #eff := 2
        self.assertEqual( model.size_of("eff"), 60 )

        model.append_effect( z2, 5, 8 )              #eff := 3
        self.assertEqual( model.size_of("eff"), 100 )

        model.append_corrstruct( iG1, 2, iA, 8, corr_eff )
        model.append_traitstruct( obs_trate_1, eff_trate_1 )
        model.append_traitstruct( obs_trate_2, eff_trate_2 )

        self.assertEqual( model.size_of("var"), 4 )
        self.assertEqual( model.size_of("cor"), 64 )
        self.assertEqual( model.size_of("cor_eff"), 2 )
        self.assertEqual( model.size_of("obs_trt"), 2 )
        self.assertEqual( model.size_of("eff_trt"), 4 )

        shapes = model.shape_of("eff")
        self.assertEqual( shapes[0][0], 5 )
        self.assertEqual( shapes[0][1], 2 )
        self.assertEqual( shapes[1][0], 5 )
        self.assertEqual( shapes[1][1], 2 )
        self.assertEqual( shapes[2][0], 5 )
        self.assertEqual( shapes[2][1], 8 )
        self.assertEqual( shapes[3][0], 5 )
        self.assertEqual( shapes[3][1], 8 )
        shapes.clear()

        shapes = model.shape_of("obs")
        for e in shapes:
            self.assertEqual( e[0], 5 )
            self.assertEqual( e[1], 1 )
        shapes.clear()

        shapes = model.shape_of("res")
        for e in shapes:
            self.assertEqual( e[0], 2 )
            self.assertEqual( e[1], 2 )
        shapes.clear()

        shapes = model.shape_of("var")
        for e in shapes:
            self.assertEqual( e[0], 2 )
            self.assertEqual( e[1], 2 )
        shapes.clear()

        shapes = model.shape_of("cor")
        for e in shapes:
            self.assertEqual( e[0], 8 )
            self.assertEqual( e[1], 8 )
        shapes.clear()

        shapes = model.shape_of("cor_eff")
        for e in shapes:
            self.assertEqual( e[0], 2 )
            self.assertEqual( e[1], 1 )
        shapes.clear()

        shapes = model.shape_of("eff_trt")
        for e in shapes:
            self.assertEqual( e[0], 2 )
            self.assertEqual( e[1], 1 )
        shapes.clear()

        shapes = model.shape_of("obs_trt")
        for e in shapes:
            self.assertEqual( e[0], 1 )
            self.assertEqual( e[1], 1 )
        shapes.clear()

        s = evogen.Pcg()
        #s.append_model( model )
        #model.print()

        model.clear_residuals()
        model.clear_observations()
        model.clear_effects()
        model.clear_corrstruct()
        model.clear_traitstruct()

        self.assertEqual( model.size_of("res"), 0 )
        self.assertEqual( model.size_of("obs"), 0 )
        self.assertEqual( model.size_of("eff"), 0 )
        self.assertEqual( model.size_of("var"), 0 )
        self.assertEqual( model.size_of("cor"), 0 )
        self.assertEqual( model.size_of("cor_eff"), 0 )
        self.assertEqual( model.size_of("obs_trt"), 0 )
        self.assertEqual( model.size_of("eff_trt"), 0 )

    def testSolver( self ):

        model = evogen.Model()
        solver = evogen.Pcg()

        #prepare data
        iR = [0.0278, -0.0102, -0.0102, 0.0371]
        y1 = [ 4.5, 2.9, 3.9, 3.5, 5.0 ]
        y2 = [ 6.8, 5.0, 6.8, 6.0, 7.5 ]
        x1 = [ 1, 0, 0, 1, 0, 1, 1, 0, 1, 0 ]
        x2 = [ 1, 0, 0, 1, 0, 1, 1, 0, 1, 0 ]
        z1 = [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 ]
        z2 = [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 ]

        sol_true = [ 4.3608, 3.3972, 0.1510, -0.0154, -0.0784, -0.0102, -0.2703, 0.2758, -0.3161, 0.2438, 6.7998, 5.8803, 0.2797, -0.0076, -0.1703, -0.0126, -0.4778, 0.5173, -0.4789, 0.3920 ]
        
        iA = [ 1.833, 0.5, 0, -0.667, 0, -1, 0, 0,
                0.5, 2, 0.5, 0, -1, -1, 0, 0,
                0, 0.5, 2, 0, -1, 0.5, 0, -1,
                -0.667, 0, 0, 1.833, 0.5, 0, -1, 0,
                0, -1, -1, 0.5, 2.5, 0, -1, 0,
                -1, -1, 0.5, 0, 0, 2.5, 0, -1,
                0, 0, 0, -1, -1, 0, 2, 0,
                0, 0, -1, 0, 0, -1, 0, 2 ]
        
        iG1 = [ 0.084, -0.0378, -0.0378, 0.042 ]

        corr_eff = [ 1, 3 ]
        eff_trate_1 = [ 0, 1 ]             #trait 1
        obs_trate_1 = 0                    #trait 1
        eff_trate_2 = [ 2, 3 ]             #trait 2
        obs_trate_2 = 1                    #trait 2

        #define the model
        model.append_residual( iR, 2)
        model.append_observation( y1, 5 )  #obs := 0
        model.append_observation( y2, 5 )  #obs := 1
        model.append_effect( x1, 5, 2 )    #eff := 0
        model.append_effect( z1, 5, 8 )    #eff := 1
        model.append_effect( x2, 5, 2 )    #eff := 2
        model.append_effect( z2, 5, 8 )    #eff := 3

        model.append_corrstruct( iG1, 2, iA, 8, corr_eff )
        
        model.append_traitstruct( obs_trate_1, eff_trate_1 )
        model.append_traitstruct( obs_trate_2, eff_trate_2 )

        solver.append_model( model )

        solver.solve()

        sol = solver.get_solution()

        for i in range(0, 20):
            self.assertAlmostEqual( sol[i], sol_true[i], 3 )

        model.clear()


    def testSolver_with_IO( self ):

        sol_true = [ 4.3608, 3.3972, 0.1510, -0.0154, -0.0784, -0.0102, -0.2703, 0.2758, -0.3161, 0.2438, 6.7998, 5.8803, 0.2797, -0.0076, -0.1703, -0.0126, -0.4778, 0.5173, -0.4789, 0.3920 ]

        model = evogen.Model()
        solver = evogen.Pcg()

        dummy_type1 = int(1)
        dummy_type2 = float(1.0)
        
        model.append_effect("tests/data/z1.dat", dummy_type1)  #eff_0
        model.append_effect("tests/data/z2.dat", dummy_type1)  #eff_1
        model.append_effect("tests/data/x1.dat", dummy_type2)  #eff_2
        model.append_effect("tests/data/x2.dat", dummy_type2)  #eff_3
        
        model.append_observation("tests/data/y1.dat")          #obs_0
        model.append_observation("tests/data/y2.dat")          #obs_1

        model.append_residual("tests/data/iR.dat")

        corr_eff = [0, 1]

        model.append_corrstruct("tests/data/iG.dat", "tests/data/ainv.dat", corr_eff);

        self.assertEqual( model.size_of("eff"), 100 )
        self.assertEqual( model.size_of("obs"), 10 )

        eff_trate_1 = [ 2, 0 ]             #trait 1
        obs_trate_1 = 0                    #trait 1

        eff_trate_2 = [ 3, 1 ]             #trait 2
        obs_trate_2 = 1                    #trait 2

        model.append_traitstruct( obs_trate_1, eff_trate_1 )
        model.append_traitstruct( obs_trate_2, eff_trate_2 )

        solver.append_model( model )

        solver.solve()

        sol = solver.get_solution()

        solver.get_solution( "py_out_solution_model_1.dat" )

        for i in range(0, 20):
            self.assertAlmostEqual( sol[i], sol_true[i], 3 )

        model.clear()

    def testSolver_model_2( self ):

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

        #----------- DATA ----------------------------
        x = [
            1, 0, 0, 1, 0,
            1, 0, 0, 0, 1,
            1, 0, 0, 0, 1,
            1, 0, 0, 1, 0,
            0, 1, 0, 1, 0,
            0, 1, 0, 0, 1,
            0, 1, 0, 0, 1,
            0, 0, 1, 0, 1,
            0, 0, 1, 1, 0,
            0, 0, 1, 0, 1]

        s = [
            1, 0, 0, 0,
            1, 0, 0, 0,
            0, 0, 1, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            1, 0, 0, 0,
            0, 0, 0, 1,
            0, 0, 0, 1,
            1, 0, 0, 0,
            0, 0, 1, 0]

        w = [
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]

        u = [
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

        iA = [
            2, 0.5, 0, 0, -1, 0.5, 0, 0, -1, 0, 0, 0, 0, 0,
            0.5, 3, 1, 0, -1, -1, 0, 0, 0.5, -1, 0, 0, -1, 0,
            0, 1, 3.5, 0, 0.5, -0.5, 0.5, -1, 0, -1, -1, 0, 0, -1,
            0, 0, 0, 1.5, 0, 0.5, -1, 0, 0, 0, 0, 0, 0, 0,
            -1, -1, 0.5, 0, 2.5, 0, 0, -1, 0, 0, 0, 0, 0, 0,
            0.5, -1, -0.5, 0.5, 0, 3.5, -1, 0, -1, 0, 0, 0, 0, -1,
            0, 0, 0.5, -1, 0, -1, 3, 0.5, 0, 0, -1, -1, 0, 0,
            0, 0, -1, 0, -1, 0, 0.5, 2.5, 0, 0, 0, -1, 0, 0,
            -1, 0.5, 0, 0, 0, -1, 0, 0, 2.5, 0, 0, 0, -1, 0,
            0, -1, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
            0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 2, 0, 0,
            0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 2, 0,
            0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 2]

        corS = [
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1]

        y = [35, 20, 25, 40, 42, 22, 35, 34, 20, 40]

        iR = [0.0029]

        iG1 = [1]

        iG2 = [0.0076, 0.0034, 0.0034, 0.0126]

        iG3 = [0.025]

        _sol = [
            14.50,
            17.878691,
            15.926014,
            20.055,
            13.198704,
            0.56389,
            -1.244389,
            1.164936,
            -0.495,
            0.640,
            -0.866,
            -1.176,
            1.93,
            -0.559,
            -1.055081,
            0.379,
            0.863317,
            -2.998,
            1.764,
            0.261565,
            -1.583161,
            0.735732,
            0.585864,
            -0.507471,
            0.841041,
            1.299317,
            -0.157915,
            0.659541,
            -0.152954,
            0.915958,
            0.442008,
            0.093056,
            0.362213,
            -1.70083,
            0.415397,
            0.824915,
            0.460519]

        #---------------------------------------------

        model = evogen.Model()
        solver = evogen.Pcg()

        #----- define the model -------

        model.append_residual(iR, 2)

        model.append_observation(y, 10) #obs := 0

        model.append_effect(x, 10, 5)  #eff := 0
        model.append_effect(s, 10, 4)  #eff := 1
        model.append_effect(w, 10, 14) #eff := 2
        model.append_effect(u, 10, 14) #eff := 3

        corr_eff_1 = [3, 2]
        corr_eff_2 = [1]

        eff_trate = [0, 3, 2, 1]
        obs_trate = 0;

        model.append_corrstruct(iG2, 2, iA, 14, corr_eff_1)
        model.append_corrstruct(iG3, 1, corS, 4, corr_eff_2)

        model.append_traitstruct(obs_trate, eff_trate)

        solver.append_model(model)

        #----------------------------------

        solver.solve()

        sol = solver.get_solution()

        solver.get_solution( "py_out_solution_model_2.dat" )

        for i in range(0, 20):
            self.assertAlmostEqual( sol[i], _sol[i], 2 )

        model.clear()


    def testSolver_model_3( self ):

        # one trait setup
        # EXAMPLE: 11.2. p.183.
        # y = Xb + Za + e

        # corr: [a], identity matrix

        #----------- DATA ----------------------------
        x = [
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1 ]

        z = [
                1.357, -0.357,  0.286,  0.286, -0.286, -1.214, -0.143, 0.071, -0.143,  1.214,
                0.357, -0.357, -0.714, -0.714, -0.286,  0.786, -0.143, 0.071, -0.143, -0.786,
                0.357,  0.643,  1.286,  0.286,  0.714, -1.214, -0.143, 0.071, -0.143,  1.214,
               -0.643, -0.357,  1.286,  0.286, -0.286, -0.214, -0.143, 0.071,  0.857,  0.214,
               -0.643,  0.643,  0.286,  1.286, -0.286, -1.214, -0.143, 0.071, -0.143,  1.214,
                0.357,  0.643, -0.714,  0.286, -0.286,  0.786, -0.143, 0.071,  0.857,  0.214,
               -0.643, -0.357,  0.286,  0.286, -0.286,  0.786, -0.143, 0.071,  0.857, -0.786,
               -0.643,  0.643,  0.286, -0.714, -0.286, -0.214, -0.143, 0.071,  0.857, -0.786  ]
        
        corZ = [
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ]

        y = [9.0, 13.4, 12.7, 15.4, 5.9, 7.7, 10.2, 4.8]

        iR = [0.0041]

        iG1 = [0.1004]

        _sol = [
                9.94,
                0.0873718,
                -0.312064,
                0.263549,
                -0.0805778,
                0.110685,
                0.139673,
                -2.3004e-07,
                -1.62874e-07,
                -0.0609667,
                -0.0158181 ]

        #---------------------------------------------

        model = evogen.Model()
        solver = evogen.Pcg()

        #----- define the model -------

        model.append_residual(iR, 1)

        model.append_observation(y, 8) #obs := 0

        model.append_effect(x, 8, 1)  #eff := 0
        model.append_effect(z, 8, 10) #eff := 1

        corr_eff = [1]

        eff_trate = [0, 1]
        obs_trate = 0

        model.append_corrstruct(iG1, 1, corZ, 10, corr_eff)

        model.append_traitstruct(obs_trate, eff_trate)

        solver.append_model(model)

        #model.print()

        #----------------------------------

        solver.solve()

        sol = solver.get_solution()
        
        solver.get_solution( "py_out_solution_model_3.dat" )

        for i in range(0, 11):
            self.assertAlmostEqual( sol[i], _sol[i], 2 )        

        model.clear()

    def testSolver_model_3_using_IO( self ):

        iR = [0.0041]

        iG1 = [0.1004]

        _sol = [
                9.94,
                0.0873718,
                -0.312064,
                0.263549,
                -0.0805778,
                0.110685,
                0.139673,
                -2.3004e-07,
                -1.62874e-07,
                -0.0609667,
                -0.0158181 ]

        
        model = evogen.Model()
        solver = evogen.Pcg()

        #----- define the model -------

        model.append_residual(iR, 1)

        model.append_observation("tests/data/model_3/y.dat") #obs := 0

        model.append_effect("tests/data/model_3/x.dat", 0)  #eff := 0
        model.append_effect("tests/data/model_3/z.dat", 0.0) #eff := 1

        corr_eff = [1]

        eff_trate = [0, 1]
        obs_trate = 0

        model.append_corrstruct(iG1, 1, "tests/data/model_3/corZ.dat", corr_eff)

        model.append_traitstruct(obs_trate, eff_trate)

        solver.append_model(model)

        #model.print()

        #----------------------------------

        solver.solve()

        sol = solver.get_solution()
        
        solver.get_solution( "py_out_solution_model_3_IO.dat" )

        for i in range(0, 11):
            self.assertAlmostEqual( sol[i], _sol[i], 2 )    

        model.clear()

        # ----------------------------------

        model.append_corrstruct(iG1, 1, "I", 10, corr_eff)
        model.append_corrstruct("tests/data/model_3/iG.dat", "I", 10, corr_eff)

        #model.print();
            
        model.clear();

    def testSolver_model_4( self ):

        iR = [0.0041]

        iG1 = [0.1004]
        
        model = evogen.Model()
        solver = evogen.Pcg()

        #----- define the model -------

        model.append_residual(iR, 1)

        model.append_observation("tests/data/model_4/obs_1.dat")             #obs := 0

        #model.append_effect("tests/data/model_4/obs_489_snp_141123.txt", 0.0)  #eff := 0
        model.append_effect("tests/data/model_4/obs_489_snp_100.txt", 0.0)  #eff := 0
        model.append_effect("tests/data/model_4/fixed_1.dat", 0.0)           #eff := 1

        corr_eff = [0]

        eff_trate = [1, 0]
        obs_trate = 0

        #model.append_corrstruct(iG1, 1, "I", 141123, corr_eff)
        model.append_corrstruct(iG1, 1, "I", 100, corr_eff)

        model.append_traitstruct(obs_trate, eff_trate)

        solver.append_model(model)

        solver.solve()
        
        solver.get_solution( "py_out_solution_model_4.dat" )

        model.clear()

        # ----------------------------------

if __name__ == '__main__':
    unittest.main()