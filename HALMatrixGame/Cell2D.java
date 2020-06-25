package HALMatrixGame;
import HAL.GridsAndAgents.AgentSQ2D;

class Cell2D extends AgentSQ2D<HALMatrixGame2D> {

    int currentState;
    int nextState;
    double fitness;

    void Init(int state) {
        currentState = state;
        nextState = state;
    }

    // calculate cell's fitness by neighborhood interactions
    void UpdateCellFitness() {
        int n = G.MapOccupiedHood(G.neighborhood,Xsq(),Ysq());
        fitness = 0.0;
        for (int i = 0; i < n; i++) {
            fitness += G.payoffs[this.currentState*G.nTypes + G.GetAgent(G.neighborhood[i]).currentState];
        }
        return;
    }

    // update state (after grid Step)
    void UpdateState() {
        this.currentState = this.nextState;
    }

    // stochastic replacement weighted by neighbors' fitness
    // assumes all payoffs are positive
    public void ReplaceProportionalToFitness() {

        int nNeighbors = G.MapOccupiedHood(G.neighborhood, Isq());
        G.fitnessPDF = new double[nNeighbors];
        double sum = 0.0;

        for (int j = 0; j < nNeighbors; j++) {
            Cell2D neighbor_cell = G.GetAgent(G.neighborhood[j]);
            G.fitnessPDF[j] = neighbor_cell.fitness;
            sum += neighbor_cell.fitness;
        }


        if (sum > 0) {
            // accounts for rare case where all neighbors have zero fitness (all same type)

            for (int i = 0; i < nNeighbors; i++) {
                G.fitnessPDF[i] /= sum;
            }

            int index = G.rn.RandomVariable(G.fitnessPDF); // index is fitness array
            this.nextState = G.GetAgent(G.neighborhood[index]).currentState;
        }

    }

    // deterministic replacement with fittest neighbor (random if equal)
    public void ReplaceWithMostFit() {

        // get new neighborhood
        int nNeighbors = G.MapOccupiedHood(G.neighborhood, Isq());

        // find fittest neighbor
        double maxFitness = -1e8;
        for (int j = 0; j < nNeighbors; j++) {
            Cell2D neighbor_cell = G.GetAgent(G.neighborhood[j]);

            if (neighbor_cell.fitness > maxFitness) {
                maxFitness = neighbor_cell.fitness;
                G.indexArray.clear();
                G.indexArray.add(G.neighborhood[j]);
            } else if (neighbor_cell.fitness == maxFitness) {
                G.indexArray.add(G.neighborhood[j]);
            }
        }

        int maxIndex = G.indexArray.get(G.rn.Int(G.indexArray.size()));
        this.nextState = G.GetAgent(maxIndex).currentState;
    }
}