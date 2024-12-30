import cv2
import numpy as np
from scipy.ndimage import label, center_of_mass
import matplotlib.pyplot as plt

# Step 1: Load the HMI image (assumed to be in grayscale)
def load_image(image_path):
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if image is None:
        raise FileNotFoundError("Image file not found.")
    return image

# Step 2: Preprocess the image
def preprocess_image(image):
    # Normalize image intensity
    normalized_image = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX)
    # Apply Gaussian blur to reduce noise
    blurred_image = cv2.GaussianBlur(normalized_image, (5, 5), 0)
    return blurred_image

# Step 3: Detect sunspots using thresholding
def detect_sunspots(image):
    # Use adaptive thresholding to identify sunspots
    _, binary_image = cv2.threshold(image, 50, 255, cv2.THRESH_BINARY_INV)
    
    # Label connected components
    labeled_array, num_features = label(binary_image)
    return labeled_array, num_features

# Step 4: Calculate sunspot coordinates
def calculate_sunspot_coordinates(labeled_array):
    coordinates = center_of_mass(labeled_array, labels=labeled_array, index=np.arange(1, labeled_array.max() + 1))
    return coordinates

# Step 5: Visualize results
def visualize_results(original_image, labeled_array, coordinates):
    plt.figure(figsize=(10, 10))
    plt.imshow(original_image, cmap='gray')
    plt.title('Detected Sunspots')
    
    for coord in coordinates:
        plt.plot(coord[1], coord[0], 'ro')  # Coordinates are (row, col)
    
    plt.show()

# Main function
def main(image_path):
    try:
        # Load and preprocess image
        image = load_image(image_path)
        preprocessed_image = preprocess_image(image)

        # Detect sunspots
        labeled_array, num_features = detect_sunspots(preprocessed_image)

        # Calculate sunspot coordinates
        coordinates = calculate_sunspot_coordinates(labeled_array)

        # Print results
        print(f"Number of sunspots detected: {num_features}")
        print("Coordinates of sunspots:", coordinates)

        # Visualize results
        visualize_results(image, labeled_array, coordinates)

    except Exception as e:
        print(f"An error occurred: {e}")

# Run the program
# Replace 'hmi_image.jpg' with the path to your HMI image file
main('hmi_image.jpg')
