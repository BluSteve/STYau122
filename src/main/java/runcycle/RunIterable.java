package runcycle;

import runcycle.structs.InputInfo;
import runcycle.structs.RunOutput;
import runcycle.structs.RunnableMolecule;

import java.util.Iterator;

public class RunIterable implements Iterable<RunOutput> {
	private final InputInfo initialInfo;
	private final RunnableMolecule[] initialRms;
	private int limit;

	public RunIterable(InputInfo initialInfo, RunnableMolecule[] initialRms) {
		this.initialInfo = initialInfo;
		this.initialRms = initialRms;
	}

	public int getLimit() {
		return limit;
	}

	public void setLimit(int limit) {
		this.limit = limit;
	}

	@Override
	public Iterator<RunOutput> iterator() {
		RunIterator res = new RunIterator(initialInfo, initialRms);
		res.setLimit(limit);
		return res;
	}
}

