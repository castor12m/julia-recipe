using Pkg
using SatelliteToolboxTle
using Dates
using DelimitedFiles
using SatelliteToolboxTle
using SatelliteToolboxSgp4

# Julian Date를 DateTime으로 변환하는 함수
function julian_to_datetime(julian_date::Float64)
    julian_day = floor(julian_date + 0.5)
    f = julian_date + 0.5 - julian_day
    if julian_day >= 2299161
        alpha = Int(floor((julian_day - 1867216.25) / 36524.25))
        a = julian_day + 1 + alpha - Int(floor(alpha / 4))
    else
        a = julian_day
    end
    b = a + 1524
    c = Int(floor((b - 122.1) / 365.25))
    d = Int(floor(365.25 * c))
    e = Int(floor((b - d) / 30.6001))

    day = b - d - Int(floor(30.6001 * e)) + f
    month = e < 14 ? e - 1 : e - 13
    year = month > 2 ? c - 4716 : c - 4715

    hours = (day - floor(day)) * 24
    minutes = (hours - floor(hours)) * 60
    seconds = (minutes - floor(minutes)) * 60

    return DateTime(year, month, Int(floor(day)), Int(floor(hours)), Int(floor(minutes)), Int(round(seconds)))
end


println("test")

tle = tle"""
MOLNIYA 2-9
1 07276U 74026A   24311.92737453  .00000129  00000+0  00000+0 0  9992
2 07276  64.3142 108.4262 6496720 269.2092  21.5970  2.45095859270797
"""
println(tle)
# >>>
# TLE: MOLNIYA 2-9 (Epoch = 2024-11-06T22:15:25.159)

# sgp4 함수 호출
sgp4d = sgp4_init(tle)
# println(sgp4d)
# >>>
# Sgp4Propagator{Float64, Float64}(2.46062142737453e6, 0.010694324306384442, 0.649672, 1.122494545786137, 1.8923941854258757, 4.698586916659919, 0.3769387585532154, 0.0, 0.0, 3.643061002965466, 0.649672, 1.122494545786137, 1.8923941854258757, 4.698586916659919, 0.3769387585532154, 0.010694974154408133, 3.643061002965466, 0.010694974154408133, 1.0, 1.880276800610929e-9, 0.7602146357549294, 0.38010792932325127, 0.8996374431623194, 0.9011844703073119, 0.43343575126762424, 0.18786655047692982, 2.53215306e-6, 0.000541314994525, 6.0412035375e-7, 0.0, 1.0353701020893047e-15, 2.0676411178658775e-9, 2.8982775111242456e-7, 0.0, 0.0, 0.0, 0.010694324299739542, -1.1969970067833648e-7, -1.6975908092549354e-6, :sdp4, Sgp4Constants{Float64}(6378.137, 0.07436685316871385, 0.00108262998905, -2.53215306e-6, -1.61098761e-6), SatelliteToolboxSgp4.Sgp4DeepSpace{Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 9.404960783147835e-8, -2.9464076367357507e-8, -2.7559260845207015e-8, 2.9294470122150765e-8, 2.3962138858665806e-9, 0.0, 0.0, 0.0, 0.3573074725804886, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.293125166389196, 0.0015368872171057342, -0.0010077886588507287, -0.0010426926308640414, -0.000494122728946911, -0.0032163094592292593, 0.0004441475295686838, 0.00023198126268007919, 0.0008624953459994163, 0.0011685390254637315, -6.400344794487877e-5, -0.001188727450896396, -0.002454759311182098, 4.859017045356936, -0.00021620676318565149, 0.00022435748745326553, 0.00016914463841741426, 6.384473773819658e-5, 0.0005899573811292837, -0.000257815993521081, 0.00012212457733399199, -0.00019126642823017668, -0.00010161508068180424, -3.369407484847565e-5, 0.00029990966186915736, 0.00038323176530455946, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, false, false, false))

# # TLE 에포크 계산 (DateTime 형식)
epoch_julian = tle_epoch(tle)
epoch_date = julian_to_datetime(epoch_julian)
# println(epoch_date)

# # 계산하려는 시간: 2024년 11월 7일 0시 0분 0초
target_date = DateTime(2024, 11, 7, 0, 0, 0)
# println(target_date)

# # 시간 차이 계산 (초 단위)
#time_difference_min = Int((Dates.DateTime(target_date) - Dates.DateTime(epoch_date)) / Millisecond(1) * (1 / 1000))
time_difference_min = (Dates.DateTime(target_date) - Dates.DateTime(epoch_date)) / Millisecond(1) * (1 / 1000) * (1 / 60)
println(time_difference_min)
# # 6275


# Afterward, we can propagate the orbit using the function sgp4!(sgp4d, t) 
# that propagates the mean elements defined in sgp4d by t minutes. 
# This function returns the position [km] and velocity [km/s] vectors 
# represented in the True Equator, Mean Equinox (TEME) reference frame.
# ex:
# # Propagate the orbit for 10 minutes.
# julia> r_teme, v_teme = sgp4!(sgp4d, 10)
r_teme, v_teme = sgp4!(sgp4d, time_difference_min)

# # r_teme: TEME(지구 중심 관성좌표계)에서의 위성 위치 벡터 (km 단위).
println("위치 (r_teme): ", r_teme)
# # v_teme: TEME 좌표계에서의 위성 속도 벡터 (km/s 단위).
println("속도 (v_teme): ", v_teme)


# for i in 1:5
#     #r_teme, v_teme = sgp4!(sgp4d, time_difference_min)

#     # 분단위로 변화해서 ㄴ넣어야함
#     r_teme, v_teme = sgp4!(sgp4d, 0.10 * (i - 1))
    
#     # # r_teme: TEME(지구 중심 관성좌표계)에서의 위성 위치 벡터 (km 단위).
#     println("위치 (r_teme): ", r_teme)
#     # # v_teme: TEME 좌표계에서의 위성 속도 벡터 (km/s 단위).
#     println("속도 (v_teme): ", v_teme)
# end
    
