#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>
#include <iomanip>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <Math/PxPyPzM4D.h>
#include <Math/Boost.h>
#include <Math/Vector4D.h>

// Function to perform two-body decay
std::pair<ROOT::Math::PxPyPzMVector, ROOT::Math::PxPyPzMVector>
two_body_decay(const ROOT::Math::PxPyPzMVector &parent, double mass1, double mass2, TRandom3 &rand_gen) {
    // Check if the decay is kinematically possible
    if (parent.M() < (mass1 + mass2)) {
        throw std::invalid_argument("Decay is kinematically forbidden.");
    }

    // Parent mass and momentum
    double m_parent = parent.M();
    double p = std::sqrt((m_parent * m_parent - (mass1 + mass2) * (mass1 + mass2)) *
                         (m_parent * m_parent - (mass1 - mass2) * (mass1 - mass2))) / (2 * m_parent);

    // Generate random decay direction (isotropic in the parent rest frame)
    double cos_theta = 2.0 * rand_gen.Rndm() - 1.0; // random cos(theta) in [-1, 1]
    double phi = 2.0 * M_PI * rand_gen.Rndm();      // random phi in [0, 2*pi]
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

    // Momentum vectors in the parent rest frame
    ROOT::Math::PxPyPzMVector particle1(
        p * sin_theta * std::cos(phi),
        p * sin_theta * std::sin(phi),
        p * cos_theta,
        mass1);

    ROOT::Math::PxPyPzMVector particle2(
        -particle1.Px(),
        -particle1.Py(),
        -particle1.Pz(),
        mass2);
    // Boost daughters back to the lab frame if the parent is moving
    ROOT::Math::Boost boost_vector(parent.BoostToCM());
    particle1 = boost_vector(particle1);
    particle2 = boost_vector(particle2);

    return {particle1, particle2};
}

// Function to smear particle momentum
ROOT::Math::PxPyPzMVector smear_momentum(const ROOT::Math::PxPyPzMVector &particle, TRandom3 &rand_gen) {
    double fluct = rand_gen.Gaus(1.0, 0.05); // Gaussian with mean=1 and sigma=0.05 (5% resolution)
    return ROOT::Math::PxPyPzMVector(
        particle.Px() * fluct,
        particle.Py() * fluct,
        particle.Pz() * fluct,
        particle.M());
}

int main(int argc, char** argv) {
    // Initialize random number generator with a seed for reproducibility
    TRandom3 rand_gen(42);

    // Set the number of events to generate
    const int nParticles = 100000;

    // Particle masses in GeV [https://pdg.lbl.gov/]
    const double mass_pi_ch = 0.13957;
    const double mass_k_zero = 0.497611;
    const double mass_d_zero = 1.86484;

    // Create a histogram to store generated invariant masses
    TH1F *hInvMass = new TH1F("hInvMass", "Invariant Mass;M(#pi^{+}#pi^{-}) [GeV];Events", 300, 0, 3);

    // Create a ROOT file to store the TTree
    TFile *fileout = new TFile("tracks.root", "RECREATE");
    TTree *tree = new TTree("tree", "Tree with tracks");

    // Define a vector to hold the tracks
    std::vector<ROOT::Math::PxPyPzMVector> tracks_vec;
    tree->Branch("tracks", &tracks_vec);

    // Loop over events and generate decays
    std::cout << "Generating " << nParticles << " events..." << std::endl;
    for (int i = 0; i < nParticles; ++i) {
        tracks_vec.clear();

        // Generate and decay a K0 particle into two charged pions
        ROOT::Math::PxPyPzMVector parent_k0;
        // Generate random momentum for parent K0
        double momentum_k0 = rand_gen.Exp(1.0); // Exponential distribution with mean=1.0
        double theta_k0 = rand_gen.Uniform(0, M_PI);
        double phi_k0 = rand_gen.Uniform(0, 2 * M_PI);
        double px_k0 = momentum_k0 * sin(theta_k0) * cos(phi_k0);
        double py_k0 = momentum_k0 * sin(theta_k0) * sin(phi_k0);
        double pz_k0 = momentum_k0 * cos(theta_k0);
        parent_k0 = ROOT::Math::PxPyPzMVector(px_k0, py_k0, pz_k0, mass_k_zero);

        // Decay K0
        try {
            auto daughters_k0 = two_body_decay(parent_k0, mass_pi_ch, mass_pi_ch, rand_gen);
            // Smear daughter momenta
            daughters_k0.first = smear_momentum(daughters_k0.first, rand_gen);
            daughters_k0.second = smear_momentum(daughters_k0.second, rand_gen);
            tracks_vec.push_back(daughters_k0.first);
            tracks_vec.push_back(daughters_k0.second);
        } catch (const std::invalid_argument &e) {
            // Skip event if decay is forbidden
            continue;
        }

        // Generate and decay a D0 particle into two charged pions
        ROOT::Math::PxPyPzMVector parent_d0;
        // Generate random momentum for parent D0
        double momentum_d0 = rand_gen.Exp(1.0); // Exponential distribution with mean=1.0
        double theta_d0 = rand_gen.Uniform(0, M_PI);
        double phi_d0 = rand_gen.Uniform(0, 2 * M_PI);
        double px_d0 = momentum_d0 * sin(theta_d0) * cos(phi_d0);
        double py_d0 = momentum_d0 * sin(theta_d0) * sin(phi_d0);
        double pz_d0 = momentum_d0 * cos(theta_d0);
        parent_d0 = ROOT::Math::PxPyPzMVector(px_d0, py_d0, pz_d0, mass_d_zero);

        // Decay D0
        try {
            auto daughters_d0 = two_body_decay(parent_d0, mass_pi_ch, mass_pi_ch, rand_gen);
            // Smear daughter momenta
            daughters_d0.first = smear_momentum(daughters_d0.first, rand_gen);
            daughters_d0.second = smear_momentum(daughters_d0.second, rand_gen);
            tracks_vec.push_back(daughters_d0.first);
            tracks_vec.push_back(daughters_d0.second);
        } catch (const std::invalid_argument &e) {
            // Skip event if decay is forbidden
            continue;
        }

        // Ensure exactly four tracks per event
        if (tracks_vec.size() != 4) {
            std::cerr << "Warning: Expected 4 tracks, got " << tracks_vec.size() << " tracks." << std::endl;
            continue;
        }

        // Fill invariant mass histogram for all unique pion pairs
        for (size_t itr1 = 0; itr1 < tracks_vec.size(); ++itr1) {
            for (size_t itr2 = itr1 + 1; itr2 < tracks_vec.size(); ++itr2) {
                double inv_mass = (tracks_vec[itr1] + tracks_vec[itr2]).M();
                hInvMass->Fill(inv_mass);
            }
        }

        // Fill the tree
        tree->Fill();

        // Optional: Print progress every 10%
        if ((i+1) % (nParticles / 10) == 0) {
            std::cout << "Progress: " << std::setw(2) << (100 * (i+1) / nParticles) << "%" << std::endl;
        }
    }

    // Write the TTree to file
    tree->Write();
    fileout->Close();

    std::cout << "Events generation completed. TTree saved to 'tracks.root'." << std::endl;

    // Fit the histogram
    std::cout << "Fitting the invariant mass histogram..." << std::endl;
    TF1 *fitFunc = new TF1("fitFunc",
                           "[0]/(sqrt(2*TMath::Pi())*[2])*TMath::Exp(-0.5*((x-[1])/[2])^2)"   // First Gaussian
                           " + [3]/(sqrt(2*TMath::Pi())*[5])*TMath::Exp(-0.5*((x-[4])/[5])^2)" // Second Gaussian
                           " + [6] + [7]*x + [8]*x^2 + [9]*x^3 + [10]*x^4 + [11]*x^5",      // Polynomial (5th degree)
                           0, 3);
    // Set initial parameters
    fitFunc->SetParameters(1, 0.5, 0.01,   // Amplitude, mean, sigma of first Gaussian
                           1, 1.85, 0.05,  // Amplitude, mean, sigma of second Gaussian
                           0, 0, 0, 0, 0);  // Polynomial coefficients

    // Set parameter limits to constrain the fit
    fitFunc->SetParLimits(1, 0.48, 0.52); // K0 mass
    fitFunc->SetParLimits(2, 0.0, 0.1);   // K0 width
    fitFunc->SetParLimits(4, 1.8, 1.9);   // D0 mass
    fitFunc->SetParLimits(5, 0.0, 0.2);   // D0 width

    // Perform the fit
    hInvMass->Fit(fitFunc, "R"); // "R" for using the range specified in TF1

    // Calculate bin width for signal event count
    double bin_width = hInvMass->GetBinWidth(1);

    // Extract and print the number of signal events for both Gaussians
    double signal1 = fitFunc->GetParameter(0) / bin_width;
    double signal1_err = fitFunc->GetParError(0) / bin_width;
    double signal2 = fitFunc->GetParameter(3) / bin_width;
    double signal2_err = fitFunc->GetParError(3) / bin_width;

    std::cout << "Number of signal events #1 (K0) = " << signal1 << " ± " << signal1_err << std::endl;
    std::cout << "Number of signal events #2 (D0) = " << signal2 << " ± " << signal2_err << std::endl;

    // Create a canvas to draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Invariant Mass", 800, 600);
    hInvMass->SetLineColor(kBlue);
    hInvMass->SetMarkerStyle(20);
    hInvMass->Draw("E");

    // Draw the fit function
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("Same");

    // Save the histogram as images
    canvas->SaveAs("invmass.pdf");
    canvas->SaveAs("invmass.png");

    std::cout << "Histogram and fit results saved as 'invmass.pdf' and 'invmass.png'." << std::endl;

    // Clean up
    delete hInvMass;
    delete fitFunc;
    delete canvas;

    return 0;
}

