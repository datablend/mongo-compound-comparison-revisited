package be.datablend.compoundcomparison.main;

import be.datablend.compoundcomparison.query.Query;
import be.datablend.compoundcomparison.setup.CompoundDatabase;

import java.io.IOException;

/**
 * User: dsuvee
 * Date: 02/02/12
 */
public class Main {

    public static void main(String[] args) throws IOException {

        // Start by creating the compound database
        CompoundDatabase compounddatabase= new CompoundDatabase();
        compounddatabase.create();

        // Execute compound similarity queries
        Query query = new Query(compounddatabase);
        // Find compounds with 10% similarity for the compound with cid 46209006
        query.findSimilarCompounds("46209006", 0.05);

    }

}
