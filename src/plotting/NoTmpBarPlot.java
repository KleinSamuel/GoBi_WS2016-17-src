package plotting;

import java.io.File;
import java.util.Vector;

import commandline.ColorCodingCommandline;
import debugStuff.DebugMessageFactory;
import io.AllroundFileWriter;
import io.ConfigReader;
import javafx.util.Pair;

public class NoTmpBarPlot extends Plot {

	Pair<Vector<Object>, Vector<Object>> pair;
	public String filename;
	public int minY = 0;
	public int maxY = 0;
	public int bottomMargin = 10;
	public int leftMargin = 5;
	public boolean logScaleY;
	private String outputPath, type;

	public NoTmpBarPlot(Pair<Vector<Object>, Vector<Object>> pair, String title, String xLab, String yLab,
			boolean logScaleY, String outputPath, String type, double maxY) {

		this.outputPath = outputPath;
		this.type = type;
		filename = title;
		this.logScaleY = logScaleY;
		setYLab(yLab);
		this.pair = pair;
		setTitle(title);
		setXLab(xLab);

		this.maxY = (int) Math.ceil(maxY);
		// maxY *= 10;
	}

	public void plot() {
		RExecutor r = new RExecutor(
				generateCommand(ConfigReader.readConfig().get("output_directory") + outputPath + "_" + type + ".png"));

		Thread t = new Thread(r);
		t.start();

		try {

			DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE,
					ColorCodingCommandline.toBlue("Wait for R to plot.."));
			t.join();
			DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE,
					"(" + this.filename + ".png) R thread terminated.");

		} catch (InterruptedException e) {
			throw new RuntimeException("R did not exit properly!");
		}
	}

	public Vector<Object> logScaleY(Vector<Object> in) {

		Vector<Object> tmp = new Vector<>();
		double max = 0;

		for (Object v : in) {
			double x = Math.log10((double) ((int) v + 1.0));
			tmp.add(x);
			max = Math.max(max, x);
		}

		this.maxY = (int) (max + 1);
		return tmp;
	}

	@Override
	String generateCommand(String filename) {

		File tmp = new File(
				ConfigReader.readConfig().get("output_directory") + "/" + this.outputPath + "_" + type + ".txt");

		AllroundFileWriter.writeVector(this.pair.getKey(), tmp);
		AllroundFileWriter.writeVector(this.pair.getValue(), tmp, true);

		String command = "";
		command += String.format("png(\'%s\',width=3.25,height=3.25,units=\'in\',res=400,pointsize=4);", filename);
		command += String.format("x<-scan(\'%s\',nlines=1,skip=0);", tmp.getAbsolutePath().replace("\\", "/"));
		command += String.format("y<-scan(\'%s\',nlines=1,skip=1,what=character());",
				tmp.getAbsolutePath().replace("\\", "/"));
		command += String.format("op<-par(mar=c(" + bottomMargin + "," + leftMargin + ",4,2)+0.1);");
		// if(logScaleY){
		command += String.format("options(scipen=10);");
		// }
		command += String.format("barplot(x,names.arg=y,col=rainbow(\'%s\'),las=2" + (logScaleY ? ",log=\'y\'" : "")
				+ ",ylim=c(" + 0 + "," + maxY + "));", this.pair.getKey().size());
		command += String.format("par(op);");
		command += String.format("title(main=\'%s\', xlab=\'%s\', ylab=\'%s\');", super.title, super.xLab, super.yLab);
		command += "dev.off();";

		return command;
	}
}