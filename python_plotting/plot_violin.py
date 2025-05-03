import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main(generation, phenotype_id):
    # Get the script and project directories
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, '..'))

    # Paths to input data
    meg_path = os.path.join(project_root, 'output', 'meg', 'meg_burst_stats_merged.csv')
    sim_path = os.path.join(project_root, 'output', 'model', f'generation_{generation}', 'stats', f'phenotype_{phenotype_id}_stats.csv')

    if not os.path.exists(meg_path):
        raise FileNotFoundError(f"MEG data not found: {meg_path}")
    if not os.path.exists(sim_path):
        raise FileNotFoundError(f"Simulation stats not found: {sim_path}")

    # Load data
    meg_df = pd.read_csv(meg_path)
    meg_df['Source'] = 'MEG'

    sim_df = pd.read_csv(sim_path)
    sim_df['Source'] = 'Simulation'

    # Combine stats
    stats_to_plot = ['burst_rate', 'mean_duration', 'mean_amplitude']
    combined = pd.concat([meg_df, sim_df], ignore_index=True)[['Source'] + stats_to_plot]

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for i, stat in enumerate(stats_to_plot):
        sns.violinplot(data=combined, x='Source', y=stat, ax=axes[i])
        axes[i].set_title(stat.replace('_', ' ').title())

    fig.suptitle(f'Burst Feature Comparison – Gen {generation}, Pheno {phenotype_id}')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Ensure output directory exists
    plot_dir = os.path.join(script_dir, 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    # Save the figure
    output_path = os.path.join(plot_dir, f'violin_bursts_gen{generation}_pheno{phenotype_id}.png')
    plt.savefig(output_path)
    plt.close()
    print(f"Saved plot to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 plot_violin.py <generation> <phenotype_id>")
        sys.exit(1)

    generation = int(sys.argv[1])
    phenotype_id = int(sys.argv[2])
    main(generation, phenotype_id)
