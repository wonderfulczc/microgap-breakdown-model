#pragma once
namespace streamer_rf::streamer {
struct TransportCoefficients {
 double mobility{},diffusion{},ionization_frequency{},attachment_two_body_frequency{},attachment_three_body_frequency{};
 double ionization_townsend{},attachment_two_body_townsend{},attachment_three_body_townsend{},drift_speed{};
};
TransportCoefficients evaluate_morrow_lowke(double electric_field,double neutral_density,double pressure,double gas_temperature);
double morrow_lowke_breakdown_field(double neutral_density,double pressure,double gas_temperature);
}
