package HALMatrixGame;
import HAL.GridsAndAgents.AgentGrid3D;
import HAL.Rand;
import static HAL.Util.*;
import HAL.Tools.FileIO;

import java.io.File;
import java.util.ArrayList;

public class HALMatrixGame3D extends AgentGrid3D<Cell3D> {

    public static final String[] labels = new String[] {"Type 1","Type 2","Type 3"};
    public static final int DETERMINISTIC = 1, STOCHASTIC = 2, nTypes = 3;
    public double UPDATE_FRACTION = 1.0;
    public int PROCESS = DETERMINISTIC;

    public static double[] payoffs = new double[]{
            0.7,0.0,0.7,
            0.3,0.4,0.8,
            1.0,0.3,0.2};

    // default names:
    public static String MESH_GRID_FILENAME = "IsomatrixGrid.csv";
    public static String SINGLE_SIMULATION_FILENAME = "HAL_trajectory.csv";
    public static String PAYOFF_FILENAME = "payoff.csv";

    // index array for tracking fitness of neighborhood (do not change)
    public ArrayList<Integer> indexArray = new ArrayList<>();

    // neighborhood (3D moore)
    int[]neighborhood = MooreHood3D(true);//MooreHood(true);
    double[] fitnessPDF; // used in stochastic updating

    Rand rn = new Rand(); // random number generator

    public HALMatrixGame3D(int sideLength) {
        super(sideLength, sideLength, sideLength, Cell3D.class, false, false, false);

        // initialize all as max type (i.e. rounds up for this type)
        for (int i = 0; i < length; i++) { NewAgentSQ(i).Init(rn.Int(3)); }
    }


    // default constructor
    public HALMatrixGame3D(int sideLength, double[] p) {
        super(sideLength, sideLength, sideLength, Cell3D.class, false, false, false);

        // normalize
        double sum = (p[0] + p[1] + p[2]);
        int maxAt = 0;

        // find largest:
        for (int i = 0; i < p.length; i++) {
            maxAt = p[i] > p[maxAt] ? i : maxAt;
            p[i] /= sum;
        }

        // initialize all as max type (i.e. rounds up for this type)
        for (int i = 0; i < length; i++) { NewAgentSQ(i).Init(maxAt); }
        ShuffleAgents(rn);

        // other two types
        int type1 = (maxAt + 1) % 3;
        int type2 = (maxAt + 2) % 3;

        int[] updateN = new int[] {(int)Math.round(p[type1]*(double)length),(int)Math.round(p[type2]*(double)length)};

        int count = 0;
        int i = 0;

        for (Cell3D c : this) {
            if (count < updateN[i]) {
                c.Init((i==0) ? type1 : type2);
                count++;
            } else {
                count = 0;
                i++;
            }

            if (i >= updateN.length) {
                break;
            }
        }
    }

    public void Output(boolean append) {
        String folder = "HALMatrixGame/HALMatrix-output/";
        File dir = new File(folder);
        dir.mkdir();

        String payoff_filename = folder + PAYOFF_FILENAME;
        String filename = folder + SINGLE_SIMULATION_FILENAME;

        FileIO fileIO= new FileIO(filename, (append) ? "a" : "w");
        StringBuilder sb = new StringBuilder();
        if (!append) {
            FileIO fileIO2= new FileIO(payoff_filename, "w");
            // output payoff matrix
            PrintPayoff(fileIO2);
            fileIO2.Close();

            if (!append) {
                // setup header:
                sb.append("time" + "," + labels[0] + "," + labels[1] + "," + labels[2] + "\n");
                fileIO.Write(sb.toString());
            }
        }

        sb = new StringBuilder();

        int[] n = new int[nTypes];
        for (int i = 0; i < length; i++) {
            n[GetAgent(i).currentState]++;
        }
        sb.append(this.GetTick() + "," + n[0]);
        for (int i = 1; i < n.length; i++) {
            sb.append("," + n[i]);
        }
        sb.append("\n");
        fileIO.Write(sb.toString());
        fileIO.Close();
    }

    public void Step(boolean output_trajectory) {

        if ((GetTick() == 0) && output_trajectory) {
            Output(false);
        }

        // update every cell's fitness:
        for (Cell3D cell : this) { cell.UpdateCellFitness(); }

        // how many cells to go through death-birth process
        int updateN = (int)Math.round(xDim*yDim * UPDATE_FRACTION);

        int count = 0;

        // sample without replacement:
        for (Cell3D focal_cell : this) {

            // one random individual chosen for death:
            if (PROCESS == STOCHASTIC) {
                focal_cell.ReplaceProportionalToFitness();
            } else {
                focal_cell.ReplaceWithMostFit();
            }

            if (count >= updateN) { break;}
            count++;
        }

        for (Cell3D c : this) { c.UpdateState(); }

        // clean, shuffle (so that next iteration is random order)
        CleanAgents();
        ShuffleAgents(rn);
        IncTick();

        if (output_trajectory) { Output(true); }

    }

    public static void PrintPayoff(FileIO file) {
        System.out.println("Payoff Matrix:");

        StringBuilder sb = new StringBuilder();
        sb.append(labels[0]+","+labels[1]+","+labels[2]+"\n");
        file.Write(sb.toString());

        for (int i = 0; i < nTypes; i++) {
            sb = new StringBuilder();

            sb.append(payoffs[i*nTypes]);

            for (int j = 1; j < nTypes; j++) {
                sb.append(", " + payoffs[i*nTypes+ j]);
            }

            // write out to screen
            System.out.println(sb);

            // write out to file
            sb.append("\n");
            file.Write(sb.toString());
        }
    }

    public void SingleSimulation(int timesteps) {

        // set up the model (with given sidelength and initial proportions)
//        HALMatrixGame3D model = new HALMatrixGame3D(sideLength,x0);

        for (int time = 0; time <= timesteps; time++) {
            System.out.println("time: " + time);
            this.Step(true);
        }
    }

    public void MeshGrid(int timesteps, int nSims) {

        UPDATE_FRACTION = 1.0;
        int sideLength = this.xDim;

        int grid_divisions=20; // mesh fine-ness

        String folder = "HALMatrixGame/HALMatrix-output/";
        File dir = new File(folder);
        dir.mkdir();

        // output payoff matrix
        String payoff_filename = folder + PAYOFF_FILENAME;
        FileIO fileIO2= new FileIO(payoff_filename, "w");
        HALMatrixGame3D.PrintPayoff(fileIO2);
        fileIO2.Close();


        String filename = folder + MESH_GRID_FILENAME;
        FileIO fileIO= new FileIO(filename,"w");
        fileIO.Close();

        double[] p = new double[] {0.0,0.0,0.0};


        for (int i = 0; i <= grid_divisions; i++) {
            for (int j = 0; j <= grid_divisions; j++) {
                for (int k = 0; k <= grid_divisions; k++) {
                    // choose only the cases where i + j + k sum to grid divisions
                    if ((i + j + k) == grid_divisions) {

                        p[0] = (i==0) ? 0.0 : (double)((double) i / (double)grid_divisions);
                        p[1] = (j==0) ? 0.0 : (double)((double) j / (double)grid_divisions);
                        p[2] = (k==0) ? 0.0 : (double)((double) k / grid_divisions);

                        System.out.println(p[0] + "," + p[1] + "," + p[2]);

                        // run model one single step and then output it:

                        for (int l = 0; l < nSims; l++) {

                            HALMatrixGame3D model = new HALMatrixGame3D(sideLength, p);
                            int[] x0 = new int[]{0, 0, 0};

                            for (Cell3D c : model) {
                                x0[c.currentState]++;
                            }

                            for (int m = 0; m < timesteps; m++) {
                                model.Step(false);
                            }

                            int[] xF = new int[]{0, 0, 0};
                            for (Cell3D c : model) {
                                xF[c.currentState]++;
                            }

                            // write out
                            fileIO = new FileIO(filename, "a");
                            StringBuilder sb = new StringBuilder();
                            sb.append(x0[0] + "," + x0[1] + "," + x0[2] + ",");
                            sb.append(xF[0] + "," + xF[1] + "," + xF[2] + "," + l + "," + grid_divisions + "\n");
                            fileIO.Write(sb.toString());
                            fileIO.Close();

                        }

                    }

                }
            }
        }
    }

}