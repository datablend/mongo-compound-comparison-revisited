package be.datablend.compoundcomparison.query;

import be.datablend.compoundcomparison.setup.CompoundDatabase;
import com.mongodb.*;

import java.io.IOException;
import java.util.*;

import static be.datablend.compoundcomparison.definition.Definition.*;

/**
 * User: dsuvee
 * Date: 02/02/12
 */
public class Query {

    private DBCollection compoundsCollection = null;
    private DBCollection fingerprintCountsCollection = null;

    public Query(CompoundDatabase compoundDatabase) {
         compoundsCollection = compoundDatabase.getCompoundsCollection();
         fingerprintCountsCollection = compoundDatabase.getFingerprintCountsCollection();
    }

    public void findSimilarCompounds(String compound, double similarity) throws IOException {

        // Retrieve the particular compound
        DBObject object = compoundsCollection.findOne(QueryBuilder.start("compound_cid").is(compound).get());

        // Retrieve the relevant properties
        String pubchemcid = (String)object.get(COMPOUNDCID_PROPERTY);
        String smiles = (String)object.get(SMILES_PROPERTY);
        System.out.println("Trying to find " + (similarity * 100) + "% similar molecues for PubChem compound " + pubchemcid  + " ( " + smiles + " )");

        // Extract the fingerprints to find
        List<Integer> fingerprintstofind = Arrays.asList(((BasicDBList)object.get(FINGERPRINTS_PROPERTY)).toArray(new Integer[] {}));

        // Sort the fingerprints on total number of occurences
        fingerprintstofind = findSortedFingerprints(fingerprintstofind);

        // Execute the query using the mongodb query functionalities
        long startnative = System.currentTimeMillis();
        executeNativeQuery(similarity, fingerprintstofind);
        long stopnative = System.currentTimeMillis();
        System.out.println("Total time for native query: " + (stopnative-startnative) + " ms\n");
        long startmr = System.currentTimeMillis();
        executeMapReduceQuery(similarity, fingerprintstofind);
        long stopmr = System.currentTimeMillis();
        System.out.println("Total time for map reduce query: " + (stopmr-startmr) + " ms\n");
        long startaggr = System.currentTimeMillis();
        executeAggregrateQuery(similarity, fingerprintstofind);
        long stopaggr = System.currentTimeMillis();
        System.out.println("Total time for map reduce query: " + (stopaggr-startaggr) + " ms\n");
    }

    private void executeNativeQuery(double similarity, List<Integer> fingerprintsToFind) {
        System.out.println("Executing using native query ...");

        // Calculate the essential numbers
        int maxnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() / similarity);
        int minnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() * similarity);
        int numberoffingerprintstoconsider = fingerprintsToFind.size() - minnumberofcompoundfingerprints;
        List<Integer> fingerprintsToConsider = fingerprintsToFind.subList(0,numberoffingerprintstoconsider+1);

        // Find all compounds that have one (or more) fingerprints in the list that we consider, but for which the total count of fingeprints is within the defined limits
        DBObject compoundquery = QueryBuilder.start(FINGERPRINTS_PROPERTY).in(fingerprintsToConsider).and(FINGERPRINTCOUNT_PROPERTY).lessThanEquals(maxnumberofcompoundfingerprints).and(FINGERPRINTCOUNT_PROPERTY).greaterThanEquals(minnumberofcompoundfingerprints).get();

        // Execute the query
        DBCursor compounds = compoundsCollection.find(compoundquery);
        int numberofresults = 0;
        // Let's process the found compounds locally
        while(compounds.hasNext()) {
            DBObject compound = compounds.next();
            // Retrieve all fingerprints
            BasicDBList fingerprints = ((BasicDBList) compound.get(FINGERPRINTS_PROPERTY));
            // Calculate the intersection on the total list of fingerprints
            fingerprints.retainAll(fingerprintsToFind);

            // If the remaining list of fingeprints contains at least the minimun number of featues
            if (fingerprints.size() >= minnumberofcompoundfingerprints) {
                // Retrieve the total count
                int totalcount = (Integer)compound.get(FINGERPRINTCOUNT_PROPERTY);
                // Calculate the tanimoto coefficient
                double tanimoto = (double) fingerprints.size() / (totalcount + fingerprintsToFind.size() - fingerprints.size());
                // Although we reduced the search space, we still need to check whether the tanimoto is really >= the required similarity
                if (tanimoto >= similarity) {
                    numberofresults++;
                }
            }
        }

        System.out.println(numberofresults + " matching compounds");
    }

    private void executeMapReduceQuery(double similarity, List<Integer> fingerprintsToFind) {
        System.out.println("Executing using map reduce query ...");

        // Calculate the essential numbers
        int maxnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() / similarity);
        int minnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() * similarity);
        int numberoffingerprintstoconsider = fingerprintsToFind.size() - minnumberofcompoundfingerprints;
        List<Integer> fingerprintsToConsider = fingerprintsToFind.subList(0,numberoffingerprintstoconsider+1);

        // Find all compounds that have one (or more) fingerprints in the list that we consider, but for which the total count of fingeprints is within the defined limits
        DBObject compoundquery = QueryBuilder.start(FINGERPRINTS_PROPERTY).in(fingerprintsToConsider).and(FINGERPRINTCOUNT_PROPERTY).lessThanEquals(maxnumberofcompoundfingerprints).and(FINGERPRINTCOUNT_PROPERTY).greaterThanEquals(minnumberofcompoundfingerprints).get();

        // The map fuction
        String map = "function() {  " +
                        "var found = 0; " +
                        "var fingerprintslength = this.fingerprints.length; " +
                        "for (i = 0; i < fingerprintslength; i++) { " +
                            "if (fingerprintstofind[this.fingerprints[i]] === true) { found++; } " +
                        "} " +
                        "if (found >= minnumberofcompoundfingerprints) { emit (this.compound_cid, {found : found, total : this.fingerprint_count, smiles: this.smiles} ); } " +
                     "}";

        // Execute the map reduce function
        MapReduceCommand mr = new MapReduceCommand(compoundsCollection, map, "", null, MapReduceCommand.OutputType.INLINE, compoundquery);

        // Create a hashmap for the fingerprints to find (to speed up the javascript execution)
        Map<Integer,Boolean> tofind = new HashMap<Integer,Boolean>();
        for(Integer fingerprinttofind : fingerprintsToFind) {
            tofind.put(fingerprinttofind,true);
        }

        // Set the map reduce scope
        Map<String,Object> scope = new HashMap<String,Object>();
        scope.put("fingerprintstofind",tofind);
        scope.put("minnumberofcompoundfingerprints",minnumberofcompoundfingerprints);
        mr.setScope(scope);

        // Execute the map reduce
        MapReduceOutput out = compoundsCollection.mapReduce(mr);
        int numberofresults = 0;
        // Iterate the results
        for (DBObject result : out.results()) {
            String compound_cid = (String)result.get("_id");
            DBObject value = (DBObject)result.get("value");

            // Calculate the tanimoto coefficient
            double totalcount = (Double)value.get("total");
            double found = (Double)value.get("found");
            double tanimoto = (Double)value.get("found") / ((Double)value.get("total") + fingerprintsToFind.size() - (Double)value.get("found"));
            // Although we reduced the search space, we still need to check whether the tanimoto is really >= the required similarity
            if (tanimoto >= similarity) {
                numberofresults++;
            }
        }

        System.out.println(numberofresults + " matching compounds");

    }

    private void executeAggregrateQuery(double similarity, List<Integer> fingerprintsToFind) {
        System.out.println("Executing using aggregation framework ...");

        // Calculate the essential numbers
        int maxnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() / similarity);
        int minnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() * similarity);
        int numberoffingerprintstoconsider = fingerprintsToFind.size() - minnumberofcompoundfingerprints;
        List<Integer> sublistToConsider = fingerprintsToFind.subList(0,numberoffingerprintstoconsider+1);

        // Retrieve the db itself to be able to execute the aggregration command
        DB db = compoundsCollection.getDB();

        // Pipeline holding the 6 pipeline operators
        DBObject[] pipeline = new DBObject[6];

        // Create the compound matcher pipe
        // Matches all compounds that have one (or more) fingerprints in the list that we consider, but for which the total count of fingeprints is within the defined limits
        DBObject matcher = BasicDBObjectBuilder.start(FINGERPRINTCOUNT_PROPERTY, BasicDBObjectBuilder.start("$gte", minnumberofcompoundfingerprints).add("$lte", maxnumberofcompoundfingerprints).get()).get();
        pipeline[0] = BasicDBObjectBuilder.start("$match",matcher).get();

        // Create the unwind pipe
        // Unwinds all fingerprints
        pipeline[1] = BasicDBObjectBuilder.start("$unwind","$fingerprints").get();

        // Create the fingerprint matcher pipe
        // Only maintain the documents from which the fingerprint is in the list of the fingerprints to find
        pipeline[2] = BasicDBObjectBuilder.start("$match",BasicDBObjectBuilder.start("fingerprints", BasicDBObjectBuilder.start("$in",fingerprintsToFind).get()).get()).get();

        // Create the fingerprint group pipe
        // Group all documents on compound id and count the number of remaining fingeprints
        DBObject grouping = BasicDBObjectBuilder.start("_id", "$compound_cid").add("fingerprintmatches", BasicDBObjectBuilder.start("$sum", 1).get()).add("totalcount", BasicDBObjectBuilder.start("$first","$fingerprint_count").get()).add("smiles", BasicDBObjectBuilder.start("$first","$smiles").get()).get();
        pipeline[3] = BasicDBObjectBuilder.start("$group",grouping).get();

        // Create the compound project pipe
        // Calculates the tanimoto coefficient
        Object[] additionarray = new Object[] { fingerprintsToFind.size(), "$totalcount" };
        DBObject addition = BasicDBObjectBuilder.start("$add", additionarray ).get();
        Object[] subtractionarray = new Object[] { addition , "$fingerprintmatches" };
        DBObject subtraction = BasicDBObjectBuilder.start("$subtract", subtractionarray ).get();
        Object[] divisionarray = new Object[] { "$fingerprintmatches", subtraction };
        DBObject tanimoto = BasicDBObjectBuilder.start("$divide", divisionarray ).get();
        DBObject projection = BasicDBObjectBuilder.start("_id",1).add("tanimoto",tanimoto).add("smiles",1).get();
        pipeline[4] = BasicDBObjectBuilder.start("$project",projection).get();

        // Create the tanimoto matcher pipe
        // Matches all compounds that satisfy the tanimoto requirement
        DBObject tanimotomatcher = BasicDBObjectBuilder.start("tanimoto", BasicDBObjectBuilder.start("$gte", similarity).get()).get();
        pipeline[5] = BasicDBObjectBuilder.start("$match",tanimotomatcher).get();

        // Create the command
        DBObject aggregatecommand = BasicDBObjectBuilder.start("aggregate", "compounds").add("pipeline",pipeline).get();
        System.out.println(aggregatecommand);
        CommandResult commandresult = db.command(aggregatecommand);

        System.out.println(((BasicDBList)commandresult.get("result")).size() + " matching compounds");

    }

    // Helper method to retrieve the fingeprints, but sorted on total count over the entire compound population
    private List<Integer> findSortedFingerprints(List<Integer> fingerprintsToFind) {
        System.out.println("Compound has " + fingerprintsToFind.size() + " unique fingerprints\n");

        List<Integer> sortedFingerprintsToFind = new ArrayList<Integer>();

        // Find all fingerprint count documents that have a fingerpint in the list of fingerprints to find
        DBObject fingerprintcountquery = QueryBuilder.start(FINGERPRINT_PROPERTY).in(fingerprintsToFind.toArray()).get();
        // Only retrieve the fingerprint string itself
        DBObject fingerprintcountselection = QueryBuilder.start(FINGERPRINT_PROPERTY).is(1).get();
        // Sort the result on count
        DBObject fingerprintcountsort = QueryBuilder.start(COUNT_PROPERTY).is(1).get();

        // Execute the query on the fingerprint counts collection
        DBCursor fingerprintcounts = fingerprintCountsCollection.find(fingerprintcountquery, fingerprintcountselection).sort(fingerprintcountsort);

        // Iterate and add them one by one (fingerprints with the smallest count will come first)
        while (fingerprintcounts.hasNext()) {
            DBObject fingerprintcount = fingerprintcounts.next();
            sortedFingerprintsToFind.add((Integer)fingerprintcount.get(FINGERPRINT_PROPERTY));
        }

        return sortedFingerprintsToFind;

    }

}
