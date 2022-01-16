package nddo;

public class Config { // serializable
	public int rhf_numIt_max = 1600;
	public int uhf_numIt_max = 3200;
	public double ediis_threshold = 1e-2, ediis_max_diff = 1e-5, diiserror_hard_limit = 1e-4;
	public double rhf_damp = 0.8, uhf_damp = 0.55;

	public double poplethiel_tolerable = 1e-5;
	public int poplethiel_batch_size = 1;
	public int rhf_thiel_limit = 1000, uhf_thiel_limit = 1000;
	public int geom_opt_hessian_interval = 7;

	public int progress_bar_interval = 30, log_increase_time = 600; // todo move to runcycle config
}
