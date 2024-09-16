from utils import *
import sys


def main():
    input_file = sys.argv[1]

    df = pd.read_csv(str(input_file), sep="\t")

    # get a index of random row from the dataframe
    random_row = df.iloc[np.random.randint(0, df.shape[0])]
    plot_numerical_distributors(metadata=df, target_row=random_row)
    plot_numerical_distributors(metadata=df)

    # plot_categorical_distributions(metadata=df, target_row=random_row)
    # plot_categorical_distributions(metadata=df)


if __name__ == "__main__":
    main()
