"""two curves varying up and down"""
import random
import os
import math
import copy
import STLgen
import bladecoordinates
import subprocess
import time
import stl
from stl import mesh
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import truncnorm
MAX_POPULATION_SIZE = 8  # maximum number of plants in a colony(or population)
# STANDARD Deviation
# SD for SHAPE
sigma_fin1 = 0.0005  # final standard deviation
sigma_ini1 = 8  # initial standard deviation
# SD for POSITION
sigma_fin2 = 0.0005  # final standard deviation
sigma_ini2 = 1.5  # initial standard deviation

Smin = 0  # min seeds produced
Smax = 3  # max seeds produced
n_mi = 3  # modulation index
iter_max = 50  # Maximum number of iterations to be done
# no of control points + [overhang,gap] + [deflection angle,fitness]
CHROMOSOME_SIZE = 4
# design space
min_dis = 0.004
max_dis = 0.010
min_prim = -16
max_prim = 5.0
min_sec = 15.0
max_sec = 33.0
# max and min lines to be read from cfd file
MIN_Line = 3
omega = 75  # angular velocity
lead1 = [-0.012054612, 0.196694492, 0]
lead2 = [0.050003115, 0.200079, 0]
trail1 = [0.036167406, 0.209908331, 0]
trail2 = [0.1, 0.200079, 0]
# filepaths
cfdstl = "/home/fmg/OpenFOAM/VAWT/cfd/cfdopenfoam/constant/triSurface/"
cu_file = "/home/fmg/OpenFOAM/VAWT/curves/"
cfd_call = "/home/fmg/OpenFOAM/VAWT/cfd/cfdopenfoam/Allrun"
gcost1 = "/home/fmg/OpenFOAM/VAWT/cfd/cfdopenfoam/simvalue/CP%i/forces1.dat"
gcost2 = "/home/fmg/OpenFOAM/VAWT/cfd/cfdopenfoam/simvalue/CP%i/forces2.dat"
gcost3 = "/home/fmg/OpenFOAM/VAWT/cfd/cfdopenfoam/simvalue/CP%i/forces3.dat"


# function to create file in another folder
def create_file(pathname, filename, openmode):
    filepath = os.path.join(pathname, filename)
    if not os.path.exists(pathname):
        os.makedirs(pathname)
    file1 = open(pathname+filename, openmode)
    return file1

# class that generates chromosomes


class Chromosome:
    def __init__(self, chrom_num=0, gen_number=0, mode=" "):
        # creating initial gene
        self._genes = np.zeros((CHROMOSOME_SIZE), dtype=float)
        self.I = (100*gen_number) + chrom_num
        if(mode == "initialise"):
            flag = False
            while (not flag):
                self._genes[0] = min_dis + (max_dis-min_dis) * random.random()
                self._genes[1] = min_prim + \
                    ((max_prim-min_prim) * random.random())
                self._genes[2] = min_sec + \
                    ((max_sec-min_sec) * random.random())
                print("SEED#:", chrom_num)
                curve = self.profile_gen()
                Chromosome.genstl(curve, cfdstl)
                time = self.cal_cost()
                print("TIME TAKEN:", time, " mins")
                temp = self.get_cost()
                if(temp == False):
                    flag = False
                else:
                    self._genes[-1] = temp
                    print("acceptable individual")
                    print("\nCOST: ", self._genes[-1])
                    #print("\nCurve Coordinates\n", curve)
                    break

    def cal_cost(self):
        start=time.time()
        subprocess.call([cfd_call,str(self.I)]) #ADD ALLRUN FILENAME HERE
        end=time.time()
        return (end-start)/60
        # handle = create_file("simvalue/CP%i/" %
        #                      self.I, "forceCoeffs.dat", "w+")
        # handle.write(str(random.randrange(1, 10)))
        # handle.close()

    def get_cost(self):
        cfd_file1 = open(gcost1 % self.I, "r")  # ADD FORCECOEFF FILE PATH HERE
        cfd_file2 = open(gcost2 % self.I, "r")
        cfd_file3 = open(gcost3 % self.I, "r")
        cost1, cost2, cost3 = 0.0, 0.0, 0.0
        list1, list2, list3 = cfd_file1.readlines(), cfd_file2.readlines(), cfd_file3.readlines()
        # actual code
        MAX_Line1, MAX_Line2, MAX_Line3 = len(list1), len(list2), len(list3)

        if(len(list1) >= MAX_Line1):
            for i in range(MIN_Line, MAX_Line1):
                cost1 += float(list[i].split(')')[4].split()[2])+float(list[i].split(')')[5].split()[2])
            for i in range(MIN_Line, MAX_Line2):
                cost2 += float(list[i].split(')')[4].split()[2])+float(list[i].split(')')[5].split()[2])
            for i in range(MIN_Line, MAX_Line3):
                cost3 += float(list[i].split(')')[4].split()[2])+float(list[i].split(')')[5].split()[2])
            cost2 = cost2/(MAX_Line2-MIN_Line)
            cost1 = cost1/(MAX_Line1-MIN_Line)
            cost3 = cost3/(MAX_Line3-MIN_Line)
            cfd_file1.close()
            return (cost1+cost2+cost3)*omega
        else:
            print("SIMULATION INCOMPLETE........REGENERATING INDIVIDUAL")
            return False

        # for debugging

        # return random.randint(1,6000)
    def get_genes(self):
        return self._genes

    def profile_gen(self):
        curve1 = bladecoordinates.blade1
        curve2 = bladecoordinates.blade2
        curve = [curve1, curve2]
        curve = Chromosome.transform_slat(curve, self._genes[0], self._genes[1], self._genes[2])  # transforming slat
        # creating stl file of slat
        Chromosome.genstl(curve, cfdstl)
        return curve

    @staticmethod
    def genstl(curve, filename):
        mesh1 = STLgen.STL_Gen(curve[0].transpose()[0], curve[0].transpose()[1], 'turbine1.stl')
        mesh2 = STLgen.STL_Gen(curve[1].transpose()[0], curve[1].transpose()[1], 'turbine2.stl')
        blade1mesh = STLgen.combine('turbine1.stl', 'turbine2.stl')
        blade2mesh = STLgen.combine('turbine1.stl', 'turbine2.stl')
        blade3mesh = STLgen.combine('turbine1.stl', 'turbine2.stl')
        # rotating blade2 and 3
        blade2mesh.rotate([0.0, 0.0, 0.5], math.radians(120))
        blade3mesh.rotate([0.0, 0.0, 0.5], math.radians(240))
        blade1mesh.save(filename+'VAWT1.stl', mode=stl.Mode.ASCII)
        blade2mesh.save(filename+'VAWT2.stl', mode=stl.Mode.ASCII)
        blade3mesh.save(filename+'VAWT3.stl', mode=stl.Mode.ASCII)
        # blademesh = mesh.Mesh(np.concatenate(
        #     [blade1mesh.data, blade2mesh.data, blade3mesh.data]))
        # blademesh.save(filename, mode=stl.Mode.ASCII)
        print("profile stl generated and saved")

    @staticmethod
    def transform_slat(target, sep_dist, prim_angle, sec_angle):
        # making blade horizontal
        target1 = copy.deepcopy(target[0].transpose())
        target2 = copy.deepcopy(target[1].transpose())
        target1[0] -= trail1[0]
        target2[0] -= lead2[0]
        target1[1] -= trail1[1]
        target2[1] -= lead2[1]
        slope1 = math.degrees(math.atan(abs((lead1[1] - trail1[1]) / (lead1[0] - trail1[0]))))
        slope2 = math.degrees(math.atan(abs((lead2[1] - trail2[1]) / (lead2[0] - trail2[0]))))
        target1 = target1.transpose()
        target2 = target2.transpose()
        target1 = STLgen.rotate(target1, slope1, [0, 0, 1])
        target2 = STLgen.rotate(target2, slope2, [0, 0, 1])
        # transforming
        target1 = STLgen.rotate(target1, prim_angle, [0, 0, 1])
        target2 = STLgen.rotate(target2, sec_angle, [0, 0, 1])
        target1 = target1.transpose()
        target2 = target2.transpose()
        # print("SSSSSSSSSSSSSSSS",target1)
        target1[1] += 0.2
        target2[1] += 0.2
        target1[0] -= sep_dist/2
        target2[0] += sep_dist/2
        target1 = target1.transpose()
        target2 = target2.transpose()
        target = [target1, target2]
        return target

    def __str__(self):
        return self._genes.__str__()
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# class that create one set of generations


class Population:
    def __init__(self, size, gen_num=0, mode=" "):
        self._chromosomes = []
        i = 0
        gen_num = str(gen_num)
        while i < size:
            self.add_chromosomes(Chromosome(i+1, int(gen_num), mode))
            i += 1

    def add_chromosomes(self, chromosome):
        self._chromosomes.append(chromosome)

    def get_chromosomes(self):
        pop_chroms_2d_array = np.zeros((len(self._chromosomes), CHROMOSOME_SIZE), dtype=float)
        #pop_chroms_2d_array = np.around(pop_chroms_2d_array, 6)
        for i in range(len(self._chromosomes)):
            pop_chroms_2d_array[i] = self._chromosomes[i].get_genes()
        return pop_chroms_2d_array
# class that helps in evolving and mutating the genes of the chromosomes


class GeneticAlgorithm:
    @staticmethod
    def reproduce(pop, iter):
        new_pop = copy.deepcopy(pop)
        worst_cost = pop._chromosomes[-1].get_genes()[CHROMOSOME_SIZE - 1]
        best_cost = pop._chromosomes[0].get_genes()[CHROMOSOME_SIZE - 1]
        sigma_iter1 = GeneticAlgorithm.std_deviation(iter, iter_max, sigma_ini1, sigma_fin1)  # for shape
        # sigma_iter2 = GeneticAlgorithm.std_deviation(iter, iter_max , sigma_ini2 , sigma_fin2) #for position
        if(best_cost != worst_cost):
            seed_num = 0
            # limiting the number of individuals that can reproduce
            for i in range(MAX_POPULATION_SIZE):
                ratio = (pop._chromosomes[i]._genes[- 1] -
                         worst_cost) / (best_cost - worst_cost)
                # number of seeds chromosome can produce on the basis of rank
                S = Smin + (Smax - Smin) * ratio
                # print("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS",int(S))
                for j in range(int(S)):
                    seed_num += 1
                    seed = Chromosome(seed_num, iter)
                    flag = False
                    while(not flag):
                        # seed._genes[0]= np.random.normal(pop._chromosomes[i].get_genes()[0], sigma_iter1)
                        # seed._genes[1]= np.random.normal(pop._chromosomes[i].get_genes()[1], sigma_iter1)
                        # seed._genes[2]= np.random.normal(pop._chromosomes[i].get_genes()[2], sigma_iter1)
                        seed._genes[0] = GeneticAlgorithm.get_truncated_normal(pop._chromosomes[i].get_genes()[0], sigma_iter1, min_dis, max_dis).rvs()
                        seed._genes[1] = GeneticAlgorithm.get_truncated_normal(pop._chromosomes[i].get_genes()[1], sigma_iter1, min_prim, max_prim).rvs()
                        seed._genes[2] = GeneticAlgorithm.get_truncated_normal(pop._chromosomes[i].get_genes()[2], sigma_iter1, min_sec, max_sec).rvs()
                        flag = True
                        print("SEED# ", seed_num)
                        curve = seed.profile_gen()
                        Chromosome.genstl(curve, cfdstl)
                        time = seed.cal_cost()
                        print("\nTIME TAKEN:", time, " mins")
                        temp = seed.get_cost()
                        if(temp == False):
                            flag = False
                        else:
                            seed._genes[-1] = temp
                            print("ACCEPTABLE INDIVIDUAL")
                            print("\nCOST: ", seed._genes[-1])
                            #print("\nCurve Coordinates\n", curve)
                            new_pop.add_chromosomes(seed)
            GeneticAlgorithm.sort(new_pop)
            for i in range(MAX_POPULATION_SIZE):
                pop._chromosomes[i] = new_pop._chromosomes[i]
        else:
            print("best and worst cost equal can`t reproduce")
            return False, False
        return pop, new_pop

    @staticmethod
    def std_deviation(iter, iter_max, sigma_ini, sigma_fin):
        sigma_iter = (((iter_max - iter)**n_mi) / iter_max **n_mi) * (sigma_ini - sigma_fin) + sigma_fin
        return sigma_iter

    @staticmethod
    def sort(pop):
        pop_chroms_2d_array = pop.get_chromosomes()
        # print("chroms",pop.get_chromosomes())
        sindices = np.argsort(pop_chroms_2d_array[:, -1], axis=0)
        print("sindices", sindices)
        sorted_chroms = Population(len(pop._chromosomes), 0, "zeros")
        for i in range(0, len(pop._chromosomes)):
            sorted_chroms._chromosomes[i]._genes = pop_chroms_2d_array[sindices[-(i+1)]]
        for i in range(0, len(pop._chromosomes)):
            pop._chromosomes[i] = sorted_chroms._chromosomes[i]

    @staticmethod
    def get_truncated_normal(mean, sd, low, upp):
        return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

#------------------------------------------------------------------------------------------------------------------------------------#-


#------------------------------------------------------------------------------------------------------------------------------------#-


def _print_population(new_pop, gen_number):
    chroms = new_pop.get_chromosomes()
    print("\n---------------------------------------------------------")
    print("PRINTING AFTER SORTING THE POPULATION")
    print("Generation#", gen_number,
          "|Fittest chromosome fitness:", chroms[0][-1])
    print("-----------------------------------------------------------")
    # dplot1.update(gen_number,chroms[0][-1][1])
    i = 0
    for i in range(MAX_POPULATION_SIZE):
        print("PLANT#", i + 1, " ", "||Fitness:", chroms[i][-1], "\n")
        print("CONTROL POINTS\n", chroms[i])
        #handle.write("PLANT NO%i")
        print("--------------------------------------------------------------")
    i = 0
    for i in range(len(new_pop._chromosomes)):
        I = (100*gen_number)+i+1
        curve_file = create_file(cu_file, "curve_%04i.dat" % I, 'w+')  # saving curve coordinates
        curve_file.write(str(new_pop._chromosomes[i]._genes)+'\n')
        curve_file.close()
        plt.figure()
        curve = new_pop._chromosomes[i].profile_gen()
        Chromosome.genstl(curve, cu_file+"cstl_%04i.stl" % I)
        plt.plot(curve[0].transpose()[0], curve[0].transpose()[1])
        plt.plot(curve[1].transpose()[0], curve[1].transpose()[1])
        plt.axes().set_aspect('equal')
        # plt.xlim(-0.04,0.05)
        # plt.ylim(-0.05,0.07)
        plt.savefig(cu_file+"fig_%04i" % I)
    print("Curve File Saved")
#-------------------------------------------------------------------------


# main
# initialising population
s = time.time()
# subprocess.call(["/home/fmg/OpenFOAM/fmg-6/run/Airfoil/ga-code/delete"])
print("Running for.......\nMAX_ITERATION:", iter_max)
print("POPULATION SIZE:", MAX_POPULATION_SIZE)
print("GENERATION#0-----------------------------------------------------------------------------------")
population = Population(MAX_POPULATION_SIZE, 0, "initialise")
# dplot1=plotter.DynamicUpdate()#plotting maximum fitness dynamically
GeneticAlgorithm.sort(population)
# handle=open(r"COST.dat",'w+')
# handle.write("Generation# 0 \n")
_print_population(population, 0)
iter = 1
sigma = [GeneticAlgorithm.std_deviation(0, iter_max, sigma_ini1, sigma_fin1)]
while iter < iter_max:
    print("**************************************EVOLUTION STARTED*********************************************************")
    sigma.append(GeneticAlgorithm.std_deviation(iter, iter_max, sigma_ini1, sigma_fin1))
    print("GENERATION#", iter,"-----------------------------------------------------------------------------------")
    print("sigma shape:", GeneticAlgorithm.std_deviation(iter, iter_max, sigma_ini1, sigma_fin1))
    print("sigma position:", GeneticAlgorithm.std_deviation(iter, iter_max, sigma_ini2, sigma_fin2))
    population, new_pop = GeneticAlgorithm.reproduce(population, iter)
    print("*************************************REPRODUCED*********************************************")
    if(population == False):
        iter += 1
        break
    _print_population(new_pop, iter)
    iter += 1
plt.figure()
plt.plot(np.arange(1, iter+1), sigma)
e = time.time()
print("TOTAL TIME TAKEN: ", ((e-s)/60), "mins")
# plt.show()
