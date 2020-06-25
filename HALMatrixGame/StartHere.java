package HALMatrixGame;

public class StartHere {

    public static final int DETERMINISTIC = 1, STOCHASTIC = 2;

    public static void main(String[] args) {

        // uncomment this line to run a single simulation example
        SingleSimulationExample();

        // uncomment this line to run the meshgrid example:
        //MeshGridExample();
    }

    // an example which can plot isomatrix diagram using "HAL_isomatrix()" MATLAB function
    public static void MeshGridExample() {
        int side_length = 20;   // domain size
        int time_steps = 1;     // how many timesteps to compute resultant velocity vector
        int nSims = 50;         // the average of how many simulations?

        // run mesh grid of simulations
        HALMatrixGame2D model = new HALMatrixGame2D(side_length);

        // for 3-dimensional games, replace the previous line with the following:
        //HALMatrixGame2D model = new HALMatrixGame2D(sideLength);

        // run the meshgrid of simulations
        model.MeshGrid(time_steps,nSims);
    }

    // an example which can plot a single trajectory isomatrix diagram using "HAL_trajectory()" MATLAB function
    public static void SingleSimulationExample() {
        // single stochastic simulation
        int sideLength = 100;
        HALMatrixGame2D matrixGame = new HALMatrixGame2D(sideLength);

        // for 3-dimensional games, replace the previous line with the following:
//        HALMatrixGame3D matrixGame = new HALMatrixGame3D(sideLength);

        // fraction of population
        matrixGame.UPDATE_FRACTION = 0.5; // must be between 0 and 1
        matrixGame.PROCESS = DETERMINISTIC; // choose DETERMINISTIC or STOCHASTIC fitness process
        matrixGame.payoffs = new double[]{
                0.7,0.0,0.7,
                0.3,0.4,0.8,
                1.0,0.3,0.2};

        // must also specify these filenames in HAL_isomatrix fncs if changed:
        matrixGame.MESH_GRID_FILENAME = "IsomatrixGrid.csv";
        matrixGame.SINGLE_SIMULATION_FILENAME = "HAL_trajectory.csv";
        matrixGame.GIF_FILENAME = "HALMatrixGame2D.gif";
        matrixGame.PAYOFF_FILENAME = "payoff.csv";

        // run simulation
        int timesteps = 100;
        matrixGame.SingleSimulation(timesteps);
    }

}