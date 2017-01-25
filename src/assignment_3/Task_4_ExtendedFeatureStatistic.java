package assignment_3;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Vector;

import bamfiles.ReadPair;
import genomeAnnotation.Chromosome;
import genomeAnnotation.Gene;
import javafx.util.Pair;
import plotting.NoTmpBarPlot;
import util.Anders;

public class Task_4_ExtendedFeatureStatistic {

	// gene matching
	double geneMatching = 0;
	// intergenic
	int intergenic = 0;
	// quality; return percentage of mapped NRP
	int mapped = 0, multimapped = 0, noMis = 0, misSmaller3 = 0, misSmaller4 = 0, noClipping = 0, clippingSmaller6 = 0;

	// intergenic; return percentage of intergenic
	int gene_proximal = 0; // geneDistance >= 500
	int antisense = 0, intergenic_spliced = 0;

	// uniqueness; return percentage of genemapping
	int gene_unique = 0, multi_gene_less_4 = 0, intergenic_antisense = 0, tr_unique = 0, gene_multi_tr = 0, merged = 0,
			merged_unique = 0, intronic = 0;

	// biotype; return
	HashMap<String, Integer> tr_unique_perBiotype;
	LinkedList<String> biotype_order;

	public Task_4_ExtendedFeatureStatistic() {
		tr_unique_perBiotype = new HashMap<>();

		HashMap<String, Integer> biotypeOcc = new HashMap<>();
		for (Chromosome c : Anders.ga.getChromosomes().values()) {
			for (Gene g : c.getAllGenesSorted()) {
				Integer i = biotypeOcc.get(g.getBiotype());
				if (i == null) {
					biotypeOcc.put(g.getBiotype(), 1);
				} else {
					biotypeOcc.put(g.getBiotype(), i + 1);
				}
			}
		}

		SortedSet<Map.Entry<String, Integer>> sortedEntries = new TreeSet<>(
				new Comparator<Map.Entry<String, Integer>>() {
					@Override
					public int compare(Map.Entry<String, Integer> e1, Map.Entry<String, Integer> e2) {
						int res = e2.getValue().compareTo(e1.getValue());
						if (res == 0) {
							res = e2.getKey().compareTo(e1.getKey());
						}
						return res;
					}
				});
		sortedEntries.addAll(biotypeOcc.entrySet());
		biotype_order = new LinkedList<>();
		int i = 0;
		Iterator<Map.Entry<String, Integer>> eIt = sortedEntries.iterator();
		Map.Entry<String, Integer> next = null;
		while (i < 10 && eIt.hasNext()) {
			next = eIt.next();
			tr_unique_perBiotype.put(next.getKey(), 0);
			i++;
			biotype_order.add(next.getKey());
		}
	}

	public void addReadPair(ReadPair rp) {
		if (!rp.getMatchedGenes().isEmpty()) {
			geneMatching++;
		}
		if (rp.isIntergenic()) {
			intergenic++;
		}
		int mis = rp.getMismatchCount();
		if (mis == 0) {
			noMis++;
		}
		if (mis < 3) {
			misSmaller3++;
		}
		if (mis < 4) {
			misSmaller4++;
		}
		if (rp.getClippingSize() == 0) {
			noClipping++;
		}
		if (rp.getClippingSize() < 6) {
			clippingSmaller6++;
		}
		if (rp.getGeneDistance() > 499) {
			gene_proximal++;
		}
		if (rp.isAntisense()) {
			antisense++;
		}
		if (rp.isIntergenic() && rp.getSplitCount() > 0) {
			intergenic_spliced++;
		}
		if (rp.getMatchedGenes().size() == 1) {
			gene_unique++;
		}
		if (rp.getMatchedGenes().size() > 1 && rp.getMatchedGenes().size() < 4) {
			multi_gene_less_4++;
		}
		if (rp.isAntisense() && rp.isIntergenic()) {
			intergenic_antisense++;
		}
		if (rp.getMatchedTranscripts() != null) {
			if (rp.getMatchedTranscripts().size() == 1) {
				tr_unique++;
				Integer count = tr_unique_perBiotype
						.get(rp.getMatchedTranscripts().getFirst().getParentalGene().getBiotype());
				if (count != null) {
					tr_unique_perBiotype.put(rp.getMatchedTranscripts().getFirst().getParentalGene().getBiotype(),
							count + 1);
				}
			}
			if (rp.getMatchedTranscripts().size() < 4 && rp.getMatchedGenes().size() == 1) {
				gene_multi_tr++;
			}
		}
		if (rp.isMerged()) {
			merged++;
			if (rp.getMatchedGenes().size() == 1) {
				merged_unique++;
			}
		}
		if (rp.isIntronic()) {
			intronic++;
		}
	}

	public void increaseMapped() {
		mapped++;
	}

	public void increaseMultimapped() {
		multimapped++;
	}

	public void plot(String mapper, String outputPath) {
		plotQuality(mapper, outputPath);
		plotIntergenic(mapper, outputPath);
		plotUniqueness(mapper, outputPath);
		plotBiotype(mapper, outputPath);
	}

	public void plotBiotype(String mapper, String outputPath) {
		Vector<Object> countVec = new Vector<>();
		for (String s : biotype_order) {
			countVec.add(((double) tr_unique_perBiotype.get(s)) / tr_unique * 100d);
		}

		Vector<Object> categorie = new Vector<>(biotype_order);

		Pair<Vector<Object>, Vector<Object>> counts = new Pair<>(countVec, categorie);

		NoTmpBarPlot p = new NoTmpBarPlot(counts, mapper + "_extendedFeatureCounts_biotypes", "categories",
				"%age of transcript-unique NRPs", false, outputPath, "extendedFeatureStatistic_biotypes", 100d);
		p.plot();
	}

	public void plotUniqueness(String mapper, String outputPath) {
		Vector<Object> countVec = new Vector<>();
		countVec.add(gene_unique / geneMatching * 100d);
		countVec.add(multi_gene_less_4 / geneMatching * 100d);
		countVec.add(intergenic_antisense / geneMatching * 100d);
		countVec.add(tr_unique / geneMatching * 100d);
		countVec.add(gene_multi_tr / geneMatching * 100d);
		countVec.add(merged / geneMatching * 100d);
		countVec.add(merged_unique / geneMatching * 100d);
		countVec.add(intronic / geneMatching * 100d);

		Vector<Object> categorie = new Vector<>();
		categorie.add("gene-unique");
		categorie.add("multi-gene");
		categorie.add("intergenic-antisense");
		categorie.add("gene-tr-unique");
		categorie.add("gene-unique-tr<4");
		categorie.add("gene-merged");
		categorie.add("gene-merged-unique");
		categorie.add("intronic");

		Pair<Vector<Object>, Vector<Object>> counts = new Pair<>(countVec, categorie);

		NoTmpBarPlot p = new NoTmpBarPlot(counts, mapper + "_extendedFeatureCounts_uniqueness", "categories",
				"%age of gene mapping NRPs", false, outputPath, "extendedFeatureStatistic_uniqueness", 100d);
		p.plot();
	}

	public void plotIntergenic(String mapper, String outputPath) {
		Vector<Object> countVec = new Vector<>();
		countVec.add(gene_proximal / ((double) intergenic) * 100);
		countVec.add(antisense / ((double) intergenic) * 100);
		countVec.add(intergenic_spliced / ((double) intergenic) * 100);

		Vector<Object> categorie = new Vector<>();
		categorie.add("gene-proximal");
		categorie.add("antisense");
		categorie.add("intergenic-spliced");

		Pair<Vector<Object>, Vector<Object>> counts = new Pair<>(countVec, categorie);

		NoTmpBarPlot p = new NoTmpBarPlot(counts, mapper + "_extendedFeatureCounts_intergenic", "categories",
				"%age of intergenic NRPs", false, outputPath, "extendedFeatureStatistic_intergenic", 100d);
		p.plot();
	}

	public void plotQuality(String mapper, String outputPath) {
		Vector<Object> countVec = new Vector<>();
		double mappedD = mapped / 2d, multimappedD = multimapped / 2d;
		countVec.add(multimappedD / mappedD * 100);
		countVec.add(noMis / mappedD * 100);
		countVec.add(misSmaller3 / mappedD * 100);
		countVec.add(misSmaller4 / mappedD * 100);
		countVec.add(noClipping / mappedD * 100);
		countVec.add(clippingSmaller6 / mappedD * 100);

		Vector<Object> categorie = new Vector<>();
		categorie.add("multimapped");
		categorie.add("no mismatch");
		categorie.add("mismatch<=2");
		categorie.add("mismatch<=3");
		categorie.add("no clipping");
		categorie.add("clipping<=5");

		Pair<Vector<Object>, Vector<Object>> counts = new Pair<>(countVec, categorie);

		NoTmpBarPlot p = new NoTmpBarPlot(counts, mapper + "_extendedFeatureCounts_quality", "categories",
				"%age of mapped NRPs", false, outputPath, "extendedFeatureStatistic_quality", 100d);
		p.plot();
	}

}
