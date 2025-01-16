import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score, median_absolute_error
import os
import argparse

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Train Ridge Regression on feature and label matrices.")
    parser.add_argument("--feature-matrix", required=True, help="Path to the feature matrix CSV file.")
    parser.add_argument("--label-matrix", required=True, help="Path to the label matrix CSV file.")
    args = parser.parse_args()

    # Paths for input and output
    data_dir = "."
    output_dir = os.path.join(data_dir, "results")
    os.makedirs(output_dir, exist_ok=True)

    # Load datasets
    labels_df = pd.read_csv(args.label_matrix)  # microRNA data
    features_df = pd.read_csv(args.feature_matrix)  # TPM-normalized gene expression data

    # Ensure no NaN or infinite values
    assert not features_df.isnull().any().any(), "Data contains NaN values in features"
    assert not labels_df.isnull().any().any(), "Data contains NaN values in labels"
    assert not np.isinf(features_df).any().any(), "Data contains infinity values in features"
    assert not np.isinf(labels_df).any().any(), "Data contains infinity values in labels"

    # Remove columns with zero variance
    features_df = features_df.loc[:, features_df.std() > 0]

    # Z-transform for features
    scaler_features = StandardScaler()
    features_df = pd.DataFrame(scaler_features.fit_transform(features_df), columns=features_df.columns)

    # Prepare dictionaries to store results
    y_data = {}
    feature_importance_data = {}
    results = []

    # Loop over all columns in labels_df (each microRNA)
    for col in labels_df.columns:
        labels = labels_df[col]

        # Split data into train and test sets (80% train, 20% test)
        X_train, X_test, y_train, y_test = train_test_split(features_df, labels, test_size=0.2, random_state=42)

        # Perform grid search without cross-validation to find the best alpha
        alphas =  [0.1, 1, 10, 100, 1000, 10000, 11000, 15000]
        best_alpha = None
        best_r2 = -np.inf

        for alpha in alphas:
            model = Ridge(alpha=alpha)
            model.fit(X_train, y_train)
            y_val_pred = model.predict(X_train)
            r2 = r2_score(y_train, y_val_pred)

            if r2 > best_r2:
                best_alpha = alpha
                best_r2 = r2

        # Train Ridge Regression model with the best alpha
        ridge_model = Ridge(alpha=best_alpha)
        ridge_model.fit(X_train, y_train)

        # Predict on test data
        y_pred = ridge_model.predict(X_test)

        # Save predictions
        y_data[f'{col}_ytest'] = y_test.tolist()
        y_data[f'{col}_ypred'] = y_pred.tolist()

        # Save feature importances
        feature_importance_data[col] = ridge_model.coef_.tolist()

        # Calculate evaluation metrics
        r2 = r2_score(y_test, y_pred)
        medae = median_absolute_error(y_test, y_pred)
        mae = mean_absolute_error(y_test, y_pred)
        y_test_mean = np.mean(y_test)
        y_pred_mean = np.mean(y_pred)
        cv_test = np.std(y_test) / y_test_mean
        cv_pred = np.std(y_pred) / y_pred_mean
        perc_25_test = np.percentile(y_test, 25)
        perc_50_test = np.percentile(y_test, 50)
        perc_75_test = np.percentile(y_test, 75)
        perc_25_pred = np.percentile(y_pred, 25)
        perc_50_pred = np.percentile(y_pred, 50)
        perc_75_pred = np.percentile(y_pred, 75)
        corr = np.corrcoef(y_test, y_pred)[0, 1]

        # Perform cross-validation (after final model is trained)
        cv_scores = cross_val_score(ridge_model, X_train, y_train, cv=5, scoring='r2', n_jobs=-1)
        cv_mean_r2 = np.mean(cv_scores)
        cv_std_r2 = np.std(cv_scores)

        # Save combined metrics
        results.append([
            col, best_alpha, r2, medae, mae, y_test_mean, y_pred_mean, cv_test, cv_pred,
            perc_25_test, perc_50_test, perc_75_test, perc_25_pred,
            perc_50_pred, perc_75_pred, corr,
            cv_mean_r2, cv_std_r2
        ])

    # Save outputs to CSV files
    pd.DataFrame(y_data).to_csv(os.path.join(output_dir, "y_test_ypred_data_combined.csv"), index=False)
    pd.DataFrame(feature_importance_data, index=features_df.columns).to_csv(os.path.join(output_dir, "weight_coefficients_all.csv"), index=True)

    columns = [
        'microRNA', 'Best Alpha', 'R-squared', 'Median Abs Err', 'Mean Abs Err', 'Mean of y_test',
        'Mean of y_pred', 'Coefficient of Variance of testing dataset',
        'Coefficient of Variance of predicted dataset', '25th percentile y_test',
        '50th percentile y_test', '75th percentile y_test', '25th percentile y_pred',
        '50th percentile y_pred', '75th percentile y_pred', 'Correlation Percentiles',
        'Mean R-squared (Cross-Validation)', 'Std R-squared (Cross-Validation)'
    ]
    pd.DataFrame(results, columns=columns).to_csv(os.path.join(output_dir, "benchmarking_results_all.csv"), index=False)

    print("Combined metrics and CV results saved.")

if __name__ == "__main__":
    main()
