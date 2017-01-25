package assignment_3;

import java.util.Vector;

import bamfiles.ReadPair;
import javafx.util.Pair;
import plotting.NoTmpBarPlot;

public class Task_2_Counts {

	private double allNRPs, mappedNRPs = 0, multimappedNRPs = 0, transcriptomicNRPs = 0, merged_trNRPs = 0,
			intronicNRPs = 0, antisenseNRPs = 0, intergenicNRPs = 0;
	private String mapper;
	private String outputPath;
	private Task_4_ExtendedFeatureStatistic extendedFeatureStats;

	public Task_2_Counts(int allNRPs, String mapper, String outputPath,
			Task_4_ExtendedFeatureStatistic extendedFeatureStats) {
		this.allNRPs = allNRPs;
		this.mapper = mapper;
		this.outputPath = outputPath;
		this.extendedFeatureStats = extendedFeatureStats;
	}

	public void addReadPair(ReadPair rp) {

		if (rp.isTranscriptomic()) {
			transcriptomicNRPs++;
		} else {
			if (rp.isMerged()) {
				merged_trNRPs++;
			} else {
				if (rp.isIntronic()) {
					intronicNRPs++;
				} else {
					if (rp.isAntisense()) {
						antisenseNRPs++;
					}
					if (rp.isIntergenic()) {
						intergenicNRPs++;
					}
				}
			}
		}
		extendedFeatureStats.addReadPair(rp);
	}

	public void increaseMultimapped() {
		multimappedNRPs++;
		extendedFeatureStats.increaseMultimapped();
	}

	public void increaseMapped() {
		mappedNRPs++;
		extendedFeatureStats.increaseMapped();
	}

	public void barPlot() {
		Vector<Object> mioReads = new Vector<>();
		mioReads.add(allNRPs / 1000000);
		mioReads.add(mappedNRPs / 2000000);
		mioReads.add(multimappedNRPs / 2000000);
		mioReads.add(transcriptomicNRPs / 1000000);
		mioReads.add(merged_trNRPs / 1000000);
		mioReads.add(intronicNRPs / 1000000);
		mioReads.add(antisenseNRPs / 1000000);
		mioReads.add(intergenicNRPs / 1000000);

		Vector<Object> categorie = new Vector<>();
		categorie.add("all\n100%");
		categorie.add("mapped\n" + String.format("%.2f%%", mappedNRPs / allNRPs * 100 / 2));
		categorie.add("multimapped\n" + String.format("%.2f%%", multimappedNRPs / allNRPs * 100 / 2));
		categorie.add("transcriptomic\n" + String.format("%.2f%%", transcriptomicNRPs / allNRPs * 100));
		categorie.add("merged_tr\n" + String.format("%.2f%%", merged_trNRPs / allNRPs * 100));
		categorie.add("intronic\n" + String.format("%.2f%%", intronicNRPs / allNRPs * 100));
		categorie.add("antisense\n" + String.format("%.2f%%", antisenseNRPs / allNRPs * 100));
		categorie.add("intergenic\n" + String.format("%.2f%%", intergenicNRPs / allNRPs * 100));

		Pair<Vector<Object>, Vector<Object>> counts = new Pair<>(mioReads, categorie);

		NoTmpBarPlot p = new NoTmpBarPlot(counts, mapper + "_featureCounts", "categories", "# NRPs in Mio", false,
				outputPath, "feature_statistic", allNRPs / 1000000);
		p.plot();
		extendedFeatureStats.plot(mapper, outputPath);
	}

	public double getAllNRPs() {
		return allNRPs;
	}

}
