package HALMatrixGame;


import static HAL.Util.CircleHood;

public class Paper {

    public static final int DETERMINISTIC = 1, STOCHASTIC = 2;

    public static void main(String[] args) {






//////////////////////////////////////////////////////////////////////////////////////////////
        // trajectories

//        // treated
//        int sideLength = 30;
//        HALMatrixGame2D matrixGame = new HALMatrixGame2D(sideLength);
//        matrixGame.payoffs = new double[]{1.2,1.0,1.0,1.1,1.0,1.025,1.4,0.7,1.1};
//
//        int type1 = 2;
//        int type2 = (type1 == 1) ? 2 : 1;
//
//        System.out.println(type1 + "," +  type2);
//
//        matrixGame.MESH_GRID_FILENAME = (type1 > type2) ? "IsomatrixGrid_1_2.csv" : "IsomatrixGrid_2_1.csv";
//        matrixGame.SINGLE_SIMULATION_FILENAME = (type1 > type2) ? "trajectory_1_2.csv" : "trajectory_2_1.csv";
//        matrixGame.GIF_FILENAME = (type1 > type2) ? "matrixGame_1_2.gif" : "matrixGame_2_1.gif";
//        matrixGame.PAYOFF_FILENAME = (type1 > type2) ? "payoff_1_2.csv" : "payoff_2_1.csv";
//
//        CircleSimulation(matrixGame, 200, type1, type2);

//////////////////////////////////////////////////////////////////////////////////////////////
        // MESHGRID

        int side_length = 50;   // domain size
        int time_steps = 3;     // how many timesteps to compute resultant velocity vector
        int nSims = 20;         // the average of how many simulations?

        HALMatrixGame2D matrixGame = new HALMatrixGame2D(side_length);
        matrixGame.payoffs = new double[]{1.2,1.0,1.0,1.1,1.0,1.025,1.4,0.7,1.1}; // treated
//        matrixGame.payoffs = new double[]{1.2,1.0,1.0,1.4,1.0,1.1,1.4,0.7,1.1}; // untreated

        matrixGame.MESH_GRID_FILENAME = "IsomatrixGrid_treated_2d.csv";
        matrixGame.PAYOFF_FILENAME = "payoff_treated_2d.csv";
        matrixGame.PROCESS = STOCHASTIC;

        // run mesh grid of simulations
        matrixGame.MeshGrid(time_steps,nSims);





    }





    public static void CircleSimulation(HALMatrixGame2D model, int timesteps, int t1, int t2) {


        // assumes constructor has already been called
        for (int i = 0; i < model.length; i++) {
            model.GetAgent(i).Init(0);
        }



        int type1 = t1;
        int type2 = t2;


        // sensitive cells
        int radius = 11;
        int[] circle = CircleHood(true, radius);
        int cells = model.MapHood(circle,model.xDim/2,model.yDim/2);
        for (int i = 0; i < cells; i++) {
            model.GetAgent(circle[i]).Init(type1);
        }

        // resistant cells
        radius = 8;
        circle = CircleHood(true, radius);
        cells = model.MapHood(circle,model.xDim/2,model.yDim/2);
        for (int i = 0; i < cells; i++) {
            model.GetAgent(circle[i]).Init(type2);
        }

        model.SingleSimulation(timesteps);

//
//        GridWindow grid = new GridWindow(sideLength,sideLength,model.SCALE_FACTOR);
//        model.Draw(grid);
//
//        for (int time = 0; time <= timesteps; time++) {
//            System.out.println("time: " + time);
//            model.Draw(grid);
//            gifMaker.AddFrame(grid);
//            model.Step(true);
//        }
//        gifMaker.Close();
//        grid.Close();
    }



}