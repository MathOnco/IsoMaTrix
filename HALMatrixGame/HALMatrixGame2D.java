package HALMatrixGame;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.Gui.GifMaker;
import HAL.Gui.GridWindow;
import HAL.Gui.UIGrid;
import HAL.Rand;
import static HAL.Util.*;
import HAL.Tools.FileIO;

import java.io.File;
import java.util.ArrayList;

public class HALMatrixGame2D extends AgentGrid2D<Cell2D> {

    public static final String[] labels = new String[] {"Strategy 1","Strategy 2","Strategy 3"};
    public static final int DETERMINISTIC = 1, STOCHASTIC = 2, nTypes = 3;

    public static double UPDATE_FRACTION = 1;
    public static int PROCESS = DETERMINISTIC; // choose DETERMINISTIC or STOCHASTIC
    public static int SCALE_FACTOR = 5; // drawing size (per lattice point)

    public static double[] payoffs = new double[]{
            0.7,0.0,0.7,
            0.3,0.4,0.8,
            1.0,0.3,0.2};

    // default names:
    public static String MESH_GRID_FILENAME = "IsomatrixGrid.csv";
    public static String SINGLE_SIMULATION_FILENAME = "HAL_trajectory.csv";
    public static String GIF_FILENAME = "HALMatrixGame2D.gif";
    public static String PAYOFF_FILENAME = "payoff.csv";

    //     must also specify in HAL_isomatrix fncs if changed:

    // index array for tracking fitness of neighborhood (do not change)
    public ArrayList<Integer> indexArray = new ArrayList<>();

    // competition neighborhood (8 nearest neighbors)
    int[]neighborhood = VonNeumannHood(true);//MooreHood(true);
    double[] fitnessPDF; // used in stochastic updating

    Rand rn = new Rand(); // random number generator

    public HALMatrixGame2D(int sideLength) {
        super(sideLength, sideLength, Cell2D.class, false, false);

        // initialize randomly:
        for (int i = 0; i < length; i++) {
            NewAgentSQ(i).Init(rn.Int(3));
        }
    }

    public HALMatrixGame2D(int sideLength, double[] p) {
        super(sideLength, sideLength, Cell2D.class, false, false);

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

        for (Cell2D c : this) {
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
        String payoff_filename = "HALMatrixGame/HALMatrix-output/" + PAYOFF_FILENAME;
        String filename = "HALMatrixGame/HALMatrix-output/" + SINGLE_SIMULATION_FILENAME;

        FileIO fileIO= new FileIO(filename, (append) ? "a" : "w");
        StringBuilder sb = new StringBuilder();
        if (!append) {
            // setup header:
            sb.append("time" + "," + labels[0] + "," + labels[1] + "," + labels[2] + "\n");
            fileIO.Write(sb.toString());
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
        for (Cell2D cell : this) { cell.UpdateCellFitness(); }

        // how many cells to go through death-birth process
        int updateN = (int)Math.round(xDim*yDim * UPDATE_FRACTION);

        int count = 0;

        // sample without replacement:
        for (Cell2D focal_cell : this) {

            // one random individual chosen for death:
            if (PROCESS == STOCHASTIC) {
                focal_cell.ReplaceProportionalToFitness();
            } else {
                focal_cell.ReplaceWithMostFit();
            }

            if (count >= updateN) { break;}
            count++;
        }

        for (Cell2D c : this) { c.UpdateState(); }

        // clean, shuffle (so that next iteration is random order)
        CleanAgents();
        ShuffleAgents(rn);
        IncTick();

        if (output_trajectory) { Output(true); }

    }

    public void Draw(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                Cell2D c = this.GetAgent(x,y);
                if (c != null) {
                    vis.SetPix(x,y,CategorialColor(c.currentState));
                }
            }
        }
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

        String folder = "HALMatrixGame/HALMatrix-output/";
        File dir = new File(folder);
        dir.mkdir();

        // output payoff matrix
        String payoff_filename = folder + PAYOFF_FILENAME;
        FileIO fileIO2= new FileIO(payoff_filename, "w");
        HALMatrixGame2D.PrintPayoff(fileIO2);
        fileIO2.Close();

        GifMaker gifMaker = new GifMaker(folder + GIF_FILENAME, 100,true);
        GridWindow grid = new GridWindow(this.xDim,this.yDim,this.SCALE_FACTOR);
        this.Draw(grid);

        for (int time = 0; time <= timesteps; time++) {
//            System.out.println("time: " + time);
            this.Draw(grid);
            gifMaker.AddFrame(grid);
            this.Step(true);
        }
        gifMaker.Close();
        grid.Close();
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
        HALMatrixGame2D.PrintPayoff(fileIO2);
        fileIO2.Close();

        String filename = folder + MESH_GRID_FILENAME;
        FileIO fileIO= new FileIO(filename,"w");

        // setup two-line header:
        fileIO = new FileIO(filename, "a");
        StringBuilder sb = new StringBuilder();
        sb.append(labels[0] + "," + labels[1] + "," + labels[2] + "," + labels[0] + "," + labels[1] + "," + labels[2] + ",-,-\n");
        sb.append("x0(1),x0(2),x0(3),xF(1),xF(2),xF(3),sim,divisions\n");
        fileIO.Write(sb.toString());


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

                            HALMatrixGame2D model = new HALMatrixGame2D(sideLength, p);
                            int[] x0 = new int[]{0, 0, 0};

                            for (Cell2D c : model) {
                                x0[c.currentState]++;
                            }

                            for (int m = 0; m < timesteps; m++) {
                                model.Step(false);
                            }

                            int[] xF = new int[]{0, 0, 0};
                            for (Cell2D c : model) {
                                xF[c.currentState]++;
                            }

                            // write out
                            fileIO = new FileIO(filename, "a");
                            sb = new StringBuilder(); // clear
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