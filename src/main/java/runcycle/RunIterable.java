package runcycle;

import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;

import java.util.Iterator;

public final class RunIterable implements Iterable<RunOutput> {
	private final RunInput runInput;
	private int limit;

	public RunIterable(RunInput runInput) {
		this.runInput = runInput;
	}

	public RunIterable(RunInput runInput, int limit) {
		this.runInput = runInput;
		this.limit = limit;
	}

	public int getLimit() {
		return limit;
	}

	public void setLimit(int limit) {
		this.limit = limit;
	}

	@Override
	public Iterator<RunOutput> iterator() {
		RunIterator res = new RunIterator(runInput);
		res.setLimit(limit);
		return res;
	}
}

