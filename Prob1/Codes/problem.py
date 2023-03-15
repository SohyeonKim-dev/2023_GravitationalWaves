import os, time, subprocess
import numpy as np
import matplotlib.pyplot as plt

class Grid:
    def __init__(self, num_cells=None, xrange=None, token=None, plot_label=None):
        if num_cells is not None and xrange is not None:
            pass
        elif token is not None:
            num_cells = int(next(token))
            xrange = (float(next(token)), float(next(token)))
        else:
            raise

        self.num_cells = num_cells
        self.xrange = xrange
        self.x, self.h = np.linspace(*xrange, num_cells + 1, retstep=True)
        self.plot_label = plot_label
    
    def __eq__(self, other):
        return self.num_cells == other.num_cells and self.xrange == other.xrange

    def __str__(self):
        return '{} {:.16e} {:.16e}'.format(
            self.num_cells,
            self.xrange[0],
            self.xrange[1]
        )

class Var:
    def __init__(self, grid, left_closed, right_closed, val=None, func=None, fill_val=None, token=None, plot_label=None):
        x = grid.x
        if not left_closed:
            x = x[1:]
        if not right_closed:
            x = x[:-1]
        num_vals = x.size

        if val is not None:
            if val.size != num_vals:
                print("Sizes are incompatible")
                print("val.size: {}".format(val.size))
                print("Var.num_vals: {}".format(num_vals))
                raise
        elif func is not None:
            val = func(x)
        elif fill_val is not None:
            val = x*0 + fill_val
        elif token is not None:
            val = np.array([float(next(token)) for _ in x])
        else:
            raise

        self.grid = grid
        self.left_closed = left_closed
        self.right_closed = right_closed
        self.x = x
        self.num_vals = num_vals
        self.val = val
        self.plot_label = plot_label

    def rms(self):
        return np.linalg.norm(self.val) / np.sqrt(self.num_vals)

    def error(self, sol_func, plot_label=None):
        return Var(self.grid, self.left_closed, self.right_closed, val=sol_func(self.x)-self.val, plot_label=plot_label)

    def labeling(self):
        plt.xlabel(self.grid.plot_label)
        plt.ylabel(self.plot_label)

    def plot(self, **kwargs):
        self.labeling()
        plt.plot(self.x, self.val, **kwargs)
    
    def scatter(self, **kwargs):
        self.labeling()
        plt.scatter(self.x, self.val, **kwargs)

    def __str__(self):
        string = str(self.grid) + '\n'
        for i, v in enumerate(self.val):
            if i != 0:
                string += '\t'
            string += "{:.16e}".format(v)
        return string

class System:
    def __init__(self, grid, sol, src, res_func):
        self.grid = grid
        self.sol = Var(grid, True, True, val=sol.val, plot_label='solution')
        self.src = Var(grid, False, False, val=src.val, plot_label='source')
        self.res = Var(grid, False, False, val=np.array([res_func(sol, src, i) for i in range(1, self.grid.num_cells)]), plot_label='residual')

    def set_error(self, sol_func):
        self.err = self.sol.error(sol_func, plot_label='error')
    
    def log10_rms_residual_h2(self):
        return np.log10(self.res.rms()*self.grid.h**2)

class Problem:
    def __init__(self, no, num_cells, xrange, timeout=10, sub_output=True, **kwargs):
        # argument check
        try:
            if num_cells < 2:
                print("num_cells should be bigger than 1")
                raise
            elif no == 3:
                if num_cells % 2 != 0:
                    print("num_cells should be even")
                    raise
            elif no == 5:
                if num_cells & (num_cells - 1) != 0:
                    print("num_cells should be power of 2")
                    raise
        except:
            print("num_cells: {}".format(num_cells))
            raise

        i = no - 1
        grid = Grid(num_cells, xrange=xrange)
        target = 'run{}'.format(no)
        capture_output = not sub_output

        self.no = no
        self.grid = grid
        self.inputfile = "input{}.txt".format(no)
        self.outputfile = "output{}.txt".format(no)

        # check source file
        src_name_c = 'main{}.c'.format(no)
        src_name_f = 'main{}.f90'.format(no)
        if os.path.exists(src_name_c):
            makefile_name = 'Makefile.c'
        elif os.path.exists(src_name_f):
            makefile_name = 'Makefile.f'
        else:
            print("There is no source file: {} or {}".format(src_name_c, src_name_f))
            raise
        
        # set Makefile
        makefile_dir = '.Makefiles/'
        if os.path.exists('main{}.c'.format(no)):
            subprocess.run(['ln', '-sf', makefile_dir + makefile_name, 'Makefile'], capture_output=True)
        else:
            print("There is no Makefiles in {}".format(makefile_dir))
            raise

        # set plot_label to grid
        grid.plot_label = [
            r'$x$',
            r'$x$',
            r'$x$',
            r'$x$',
            r'$r$'
        ][i]

        # export input file
        [
            self.export_input1,
            self.export_input2,
            self.export_input3,
            self.export_input4,
            self.export_input5
        ][i](**kwargs)

        # compile
        subprocess.run(['make', 'clean'], capture_output=True)
        try:
            subprocess.run(['make', target], check=True, capture_output=capture_output)
        except subprocess.CalledProcessError as e:
            print("[Compile Error]")
            raise

        # run
        try:
            start = time.time()
            subprocess.run('./{}'.format(target), timeout=timeout, check=True, capture_output=capture_output)
            self.time = time.time() - start
        except subprocess.TimeoutExpired:
            print("[Timeout ({}s)] ".format(timeout))
            raise
        except subprocess.CalledProcessError:
            print("[Runtime Error]")
            raise

        # import output file
        [
            self.import_output1,
            self.import_output2,
            self.import_output3,
            self.import_output4,
            self.import_output5
        ][i](**kwargs)

    def export_input1(self, src_func, bvs, **kwargs):
        self.src = Var(self.grid, False, False, func=src_func)
        with open(self.inputfile, "w") as file:
            file.write('{} {}\n'.format(*bvs))
            file.write(str(self.src))
    
    def export_input2(self, func, **kwargs):
        self.f = Var(self.grid, True, True, func=func)
        with open(self.inputfile, "w") as file:
            file.write(str(self.f))

    def export_input3(self, func, **kwargs):
        self.export_input2(func)

    def export_input4(self, src_func, bvs, **kwargs):
        self.export_input1(src_func, bvs)

    def export_input5(self, **kwargs):
        with open(self.inputfile, "w") as file:
            file.write(str(self.grid))

    def import_output1(self, bvs=None, sol_func=None, **kwargs):
        # import output
        with open(self.outputfile, "r") as file:
            token = iter(file.read().split())
        grid2 = Grid(token=token)
        sol = Var(self.grid, True, True, token=token)

        # grid check
        if self.grid != grid2:
            print("Grid info is invalid")
            print('Grid1: ' + str(self.grid))
            print('Grid2: ' + str(grid2))
            raise
        
        # boundary value check
        if self.no != 5:
            if (bvs[0], bvs[1]) != (sol.val[0], sol.val[-1]):
                print("boundary value error")
                print("1: {}".format(bvs))
                print("2: {}".format((sol.val[0], sol.val[-1])))
                raise    
        else:
            bvc1 = (4. * sol.val[1] - sol.val[2]) / (3. - self.grid.h / self.grid.xrange[0]) - sol.val[0] 
            bvc2 =  1. + self.grid.xrange[0] / self.grid.xrange[1] - sol.val[-1]
            if bvc1 != 0. or bvc2 != 0.:
                print("boundary value error")
                print("left boundary check: {}".format(bvc1))
                print("right boundary check: {}".format(bvc2))
                raise

        # residual check
        if self.no != 5:
            def res_func(sol, src, i):
                return src.val[i - 1] - (sol.val[i + 1] + sol.val[i - 1] - 2 * sol.val[i]) / sol.grid.h**2
        else:
            def res_func(sol, src, i):
                h = sol.grid.h
                x = sol.x[i]
                return - (
                    + sol.val[i - 1] - 2. * sol.val[i] + sol.val[i + 1]
                    + (- sol.val[i - 1] + sol.val[i + 1]) * h / x
                ) / h**2
            self.src = Var(self.grid, False, False, val=(self.grid.x*0.)[1:][:-1])
        system = System(self.grid, sol, self.src, res_func)
        val = system.log10_rms_residual_h2()
        if  val > - 15:
            print("residual criteria error: {}".format(val))
            raise
        
        # set error
        if self.no != 5 and sol_func is not None:
            system.set_error(sol_func)
            self.err = system.err

        self.sol = sol
        self.res = system.res
        
    def import_output2(self, func, **kwargs):
        # import output
        with open(self.outputfile, "r") as file:
            token = iter(file.read().split())
        grid = Grid(token=token)
        f2 = Var(grid, True, True, token=token)

        # grid check
        if self.grid.xrange != grid.xrange and (
            (self.no == 2 and self.grid.num_cells * 2 != grid.num_cells)
            or
            (self.no == 3 and self.grid.num_cells != grid.num_cells * 2)
        ):
                print("Grid info is invalid")
                print('Grid1: ' + str(self.grid))
                print('Grid2: ' + str(grid))
                raise

        # error of f2
        err = f2.error(func, plot_label='error')

        self.f2 = f2
        self.err = err

    def import_output3(self, func, **kwargs):
        self.import_output2(func)

    def import_output4(self, bvs=None, sol_func=None, **kwargs):
        self.import_output1(bvs=bvs, sol_func=sol_func)

    def import_output5(self, **kwargs):
        self.import_output1()