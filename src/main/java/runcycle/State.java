package runcycle;

public class State {
	private static INDDOConverter converter;

	public static INDDOConverter getConverter() {
		if (converter == null) throw new IllegalStateException("NDDOConverter not set!");

		return converter;
	}

	public static void setConverter(INDDOConverter converter) {
		if (converter == null) throw new NullPointerException("Converter cannot be null!");

		State.converter = converter;
	}
}
