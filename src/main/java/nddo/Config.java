package nddo;

public class Config { // serializable
	public int rhf_numIt_max = 1000;
	public int uhf_numIt_max = 2000;
	public double ediis_threshold = 1e-2, ediis_max_diff = 1e-5;

	public double poplethiel_ideal = 1e-10;
	public int poplethiel_batch_size = 1;
	public int rhf_thiel_limit = 1000, uhf_thiel_limit = 1000;
	public int geom_opt_hessian_interval = 7;
}
