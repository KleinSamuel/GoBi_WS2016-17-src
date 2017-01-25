package assignment_3.splicing_evidence;

import java.util.TreeMap;
import java.util.Vector;

import javafx.util.Pair;
import plotting.NoTmpLinePlot;

public class Splicing_Counts {

	private TreeMap<Integer, Integer> annotUniqZero, notAnnotUniqZero, annotUniqNotZero, notAnnotUniqNotZero;
	private TreeMap<Integer, Integer> annotNotUniqZero, notAnnotNotUniqZero, annotNotUniqNotZero,
			notAnnotNotUniqNotZero;
	private int maxX = 0, maxY = 0;

	public Splicing_Counts() {
		annotUniqZero = new TreeMap<>();
		notAnnotUniqZero = new TreeMap<>();
		annotUniqNotZero = new TreeMap<>();
		notAnnotUniqNotZero = new TreeMap<>();
		annotNotUniqZero = new TreeMap<>();
		notAnnotNotUniqZero = new TreeMap<>();
		annotNotUniqNotZero = new TreeMap<>();
		notAnnotNotUniqNotZero = new TreeMap<>();
	}

	public void addSplicingEvidence(SplicingEvidence se) {
		if (se.getAnnotated()) {
			updateCount(annotUniqZero, se.getUcount_0());
			updateCount(annotUniqNotZero, se.getUcount());
			updateCount(annotNotUniqZero, se.getCount_0());
			updateCount(annotNotUniqNotZero, se.getCount());
		} else {
			updateCount(notAnnotUniqZero, se.getUcount_0());
			updateCount(notAnnotUniqNotZero, se.getUcount());
			updateCount(notAnnotNotUniqZero, se.getCount_0());
			updateCount(notAnnotNotUniqNotZero, se.getCount());
		}

	}

	public void updateCount(TreeMap<Integer, Integer> counts, int key) {
		Integer numberOfIntrons = counts.get(key);
		if (numberOfIntrons == null) {
			counts.put(key, 1);
		} else {
			counts.put(key, numberOfIntrons + 1);
		}
	}

	public void plot(String mapper, String outputPath) {
		Vector<Vector<Object>> x = new Vector<>(), y = new Vector<>();
		x.add(mapToCumulativeVector(annotUniqZero, false));
		x.add(mapToCumulativeVector(annotUniqNotZero, false));
		x.add(mapToCumulativeVector(annotNotUniqZero, false));
		x.add(mapToCumulativeVector(annotNotUniqNotZero, false));
		x.add(mapToCumulativeVector(notAnnotUniqZero, false));
		x.add(mapToCumulativeVector(notAnnotUniqNotZero, false));
		x.add(mapToCumulativeVector(notAnnotNotUniqZero, false));
		x.add(mapToCumulativeVector(notAnnotNotUniqNotZero, false));

		y.add(mapToCumulativeVector(annotUniqZero, true));
		y.add(mapToCumulativeVector(annotUniqNotZero, true));
		y.add(mapToCumulativeVector(annotNotUniqZero, true));
		y.add(mapToCumulativeVector(annotNotUniqNotZero, true));
		y.add(mapToCumulativeVector(notAnnotUniqZero, true));
		y.add(mapToCumulativeVector(notAnnotUniqNotZero, true));
		y.add(mapToCumulativeVector(notAnnotNotUniqZero, true));
		y.add(mapToCumulativeVector(notAnnotNotUniqNotZero, true));
		Pair<Vector<Vector<Object>>, Vector<Vector<Object>>> pair = new Pair<>(x, y);

		Vector<Object> legends = new Vector<>();
		legends.add("annotUniq0");
		legends.add("annotUniq");
		legends.add("annot0");
		legends.add("annot");
		legends.add("notAnnotUniq0");
		legends.add("notAnnotUniq");
		legends.add("notAnnot0");
		legends.add("notAnnot");

		NoTmpLinePlot lp = new NoTmpLinePlot(pair, mapper + "_splicingCounts", "# counts", "# introns", maxX, maxY,
				true, true, outputPath, "splicingCounts");
		lp.addLegendVector(legends);
		lp.plot();

		try {
			Thread.sleep(1000);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	public Vector<Object> mapToCumulativeVector(TreeMap<Integer, Integer> in, boolean value) {
		Vector<Object> ret = new Vector<>();
		if (value) {
			int sum = 0;
			for (Integer i : in.values()) {
				sum += i;
				ret.add(sum);
				if (sum > maxY) {
					maxY = sum;
				}
			}
		} else {
			for (Integer i : in.keySet()) {
				ret.add(i);
			}
			if (in.lastKey() > maxX) {
				maxX = in.lastKey();
			}
		}
		return ret;
	}

}
