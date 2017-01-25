package assignment_3;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.Vector;

import bamfiles.ReadPair;
import genomeAnnotation.Gene;
import genomeAnnotation.Transcript;
import io.ConfigReader;
import javafx.util.Pair;
import plotting.NoTmpLinePlot;
import util.Anders;
import util.FPKMUtils;

public class Task_3_Gene_Counts_Container {

	// geneId, Gene_counts
	private TreeMap<Gene, Task_3_Gene_Counts> countsPerGene;
	private String ebv_mapper;
	private Task_2_Counts overallCounts;
	private String mapper;

	public Task_3_Gene_Counts_Container(String ebv_mapper, Task_2_Counts overallCounts) {
		countsPerGene = new TreeMap<>();
		this.ebv_mapper = ebv_mapper;
		this.overallCounts = overallCounts;
		mapper = ebv_mapper.substring(ebv_mapper.lastIndexOf("/") + 1);
	}

	public void addReadPair(ReadPair rp) {
		Task_3_Gene_Counts c;
		for (Gene g : rp.getSpanningGenes()) {
			c = countsPerGene.get(g);
			if (c == null) {
				c = new Task_3_Gene_Counts();
				countsPerGene.put(g, c);
			}
			c.increaseNRPwithinGenReg(rp.getPcrIndex() == 0);
			if (rp.getSpanningGenes().size() == 1) {
				c.increaseNRPwithinGenRegWithout(rp.getPcrIndex() == 0);
			}
		}
		for (Gene g : rp.getMatchedGenes()) {
			c = countsPerGene.get(g);
			c.addReadPair(rp);
			if (rp.getMatchedTranscripts() != null) {
				for (Transcript t : rp.getMatchedTranscripts()) {
					if (t.getParentalGene().equals(g)) {
						c.addTranscript(t);
						if (rp.getPcrIndex() == 0) {
							c.addNURTranscript(t);
						}
					}
				}
			}
		}
	}

	public void writeGeneCountsTSV() {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(
					new File(ConfigReader.readConfig().get("output_directory") + "/" + ebv_mapper + "_counts.tsv")));
			bw.write(
					"geneId\tNRPmapping\tNRPmappingUnique\tNRPintronic\tNRPtranscriptomic\tNRPmerged-transcriptomic\tNRPhighest_transcript\tNURmapping\tNURmappingUnique\tNURintronic\tNURtranscriptomic\tNURmerged-transcriptomic\tNURhighest_transcript\n");
			for (Entry<Gene, Task_3_Gene_Counts> e : countsPerGene.entrySet()) {
				bw.write(e.getKey().getId() + "\t" + e.getValue().toString() + "\n");
			}
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private TreeMap<Double, Integer> trFpkmOcc, intronFpkmOcc, mainTrOcc;

	public void calculateStats() {
		Anders.ga.iterator()
		trFpkmOcc = new TreeMap<>();
		intronFpkmOcc = new TreeMap<>();
		mainTrOcc = new TreeMap<>();

		double trFpkm, intronFpkm, mainTr;
		Integer tr = null, intron = null, main = null;
		for (Entry<Gene, Task_3_Gene_Counts> e : countsPerGene.entrySet()) {

			trFpkm = FPKMUtils.tr_fpkm(e.getKey(), e.getValue(), overallCounts.getAllNRPs());
			intronFpkm = FPKMUtils.intron_fpkm(e.getKey(), e.getValue(), overallCounts.getAllNRPs());
			mainTr = FPKMUtils.covP(e.getValue());

			tr = trFpkmOcc.get(trFpkm);
			if (tr == null) {
				trFpkmOcc.put(trFpkm, 1);
			} else {
				trFpkmOcc.put(trFpkm, tr + 1);
			}

			intron = intronFpkmOcc.get(intronFpkm);
			if (intron == null) {
				intronFpkmOcc.put(intronFpkm, 1);
			} else {
				intronFpkmOcc.put(intronFpkm, intron + 1);
			}

			main = mainTrOcc.get(mainTr);
			if (main == null) {
				mainTrOcc.put(mainTr, 1);
			} else {
				mainTrOcc.put(mainTr, main + 1);
			}
		}
	}

	public void plot() {
		calculateStats();
		plotTr_FPKM();
		plotIntron_FPKM();
		plotCovP();
	}

	public void plotTr_FPKM() {
		Vector<Object> x = mapToCumulativeVector(trFpkmOcc, false);
		Vector<Object> y = mapToCumulativeVector(trFpkmOcc, true);
		Vector<Vector<Object>> xs = new Vector<>();
		xs.add(x);
		Vector<Vector<Object>> ys = new Vector<>();
		ys.add(y);
		Pair<Vector<Vector<Object>>, Vector<Vector<Object>>> p = new Pair<>(xs, ys);
		NoTmpLinePlot l = new NoTmpLinePlot(p, "transcript_fpkm_" + mapper, "tr_fpkm", "# occurences",
				trFpkmOcc.lastKey(), (Integer) y.lastElement(), false, false, ebv_mapper, "tr_fpkm");
		l.plot();
	}

	public void plotIntron_FPKM() {
		Vector<Object> x = mapToCumulativeVector(intronFpkmOcc, false);
		Vector<Object> y = mapToCumulativeVector(intronFpkmOcc, true);
		Vector<Vector<Object>> xs = new Vector<>();
		xs.add(x);
		Vector<Vector<Object>> ys = new Vector<>();
		ys.add(y);
		Pair<Vector<Vector<Object>>, Vector<Vector<Object>>> p = new Pair<>(xs, ys);
		NoTmpLinePlot l = new NoTmpLinePlot(p, "intron_fpkm_" + mapper, "intron_fpkm", "# occurences",
				intronFpkmOcc.lastKey(), (Integer) y.lastElement(), false, false, ebv_mapper, "intron_fpkm");
		l.plot();
	}

	public void plotCovP() {
		Vector<Object> x = mapToCumulativeVector(mainTrOcc, false);
		Vector<Object> y = mapToCumulativeVector(mainTrOcc, true);
		Vector<Vector<Object>> xs = new Vector<>();
		xs.add(x);
		Vector<Vector<Object>> ys = new Vector<>();
		ys.add(y);
		Pair<Vector<Vector<Object>>, Vector<Vector<Object>>> p = new Pair<>(xs, ys);
		NoTmpLinePlot l = new NoTmpLinePlot(p, "main_tr_" + mapper, "main_tr_coverage", "# occurences",
				mainTrOcc.lastKey(), (Integer) y.lastElement(), false, false, ebv_mapper, "main_tr");
		l.plot();
	}

	public Vector<Object> mapToCumulativeVector(TreeMap<Double, Integer> in, boolean value) {
		Vector<Object> ret = new Vector<>();
		if (value) {
			int sum = 0;
			for (Integer i : in.values()) {
				sum += i;
				ret.add(sum);
			}
		} else {
			for (Double d : in.keySet()) {
				ret.add(d);
			}
		}
		return ret;
	}

}
