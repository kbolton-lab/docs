package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link DetermineGermlineContigPloidy}.
 *
 * The test runs the CLI tool in Cohort and Case run-modes on a small simulated data.
 *
 */
public final class DetermineGermlineContigPloidyIntegrationTest extends CommandLineProgramTest {
    private static final String GCNV_SIM_DATA_DIR = toolsTestDir + "copynumber/gcnv-sim-data/";
    private static final File TEST_CONTIG_PLOIDY_PRIOR_FILE =
            new File(GCNV_SIM_DATA_DIR, "contig_ploidy_prior.tsv");
    private static final File[] TEST_COUNT_FILES = IntStream.range(0, 20)
            .mapToObj(n -> new File(GCNV_SIM_DATA_DIR, String.format("SAMPLE_%03d_counts.tsv", n)))
            .toArray(File[]::new);
    private static final File OUTPUT_DIR = createTempDir("test-ploidy");

    @Test(groups = {"python"})
    public void testCohort() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES).forEach(argsBuilder::addInput);
        argsBuilder.add(DetermineGermlineContigPloidy.CONTIG_PLOIDY_PRIORS_FILE_LONG_NAME,
                TEST_CONTIG_PLOIDY_PRIOR_FILE)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-ploidy-cohort");
        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"}, expectedExceptions = UserException.BadInput.class)
    public void testCohortWithoutContigPloidyPriors() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES).forEach(argsBuilder::addInput);
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-ploidy-cohort");
        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"}, expectedExceptions = IllegalArgumentException.class)
    public void testCohortWithSingleSample() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addInput(TEST_COUNT_FILES[0]);
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-ploidy-cohort");
        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"}, expectedExceptions = IllegalArgumentException.class)
    public void testCohortDuplicateFiles() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES).forEach(argsBuilder::addInput);
        argsBuilder.addInput(TEST_COUNT_FILES[0]);  //duplicate
        argsBuilder.add(DetermineGermlineContigPloidy.CONTIG_PLOIDY_PRIORS_FILE_LONG_NAME,
                TEST_CONTIG_PLOIDY_PRIOR_FILE)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-ploidy-cohort");
        runCommandLine(argsBuilder);
    }

    /**
     * Use the first 5 samples as case and use the contig-ploidy model generated by {@link #testCohort()}
     */
    @Test(groups = {"python"}, dependsOnMethods = "testCohort")
    public void testCase() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES, 0, 5).forEach(argsBuilder::addInput);
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-ploidy-case")
                .add(CopyNumberStandardArgument.MODEL_LONG_NAME,
                        new File(OUTPUT_DIR, "test-ploidy-cohort-model").getAbsolutePath());
        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"}, dependsOnMethods = "testCohort", expectedExceptions = UserException.BadInput.class)
    public void testCaseWithContigPloidyPrior() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES, 0, 5).forEach(argsBuilder::addInput);
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(DetermineGermlineContigPloidy.CONTIG_PLOIDY_PRIORS_FILE_LONG_NAME,
                        TEST_CONTIG_PLOIDY_PRIOR_FILE)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-ploidy-case")
                .add(CopyNumberStandardArgument.MODEL_LONG_NAME,
                        new File(OUTPUT_DIR, "test-ploidy-cohort-model").getAbsolutePath());
        runCommandLine(argsBuilder);
    }
}
