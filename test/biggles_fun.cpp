#include "biggles/detail/fun.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

void test_proc1() {
    diag("count");
    diag("%s", count2str(1).c_str());
    diag("%s", count2str(12).c_str());
    diag("%s", count2str(123).c_str());
    diag("%s", count2str(1234).c_str());
    diag("%s", count2str(12345).c_str());
    diag("%s", count2str(123456).c_str());
    diag("%s", count2str(1234567).c_str());
    diag("%s", count2str(12345678).c_str());
    diag("%s", count2str(123456789).c_str());
    diag("negative float");
    diag("%s", float2str(-1.f).c_str());
    diag("%s", float2str(-12.f).c_str());
    diag("%s", float2str(-123.f).c_str());
    diag("%s", float2str(-1234.f).c_str());
    diag("%s", float2str(-12345.f).c_str());
    diag("%s", float2str(-123456.f).c_str());
    diag("%s", float2str(-1234567.f).c_str());
    diag("%s", float2str(-12345678.f).c_str());
    diag("%s", float2str(-123456789.f).c_str());
    diag("positive float");
    diag("%s", float2str(1.f).c_str());
    diag("%s", float2str(12.f).c_str());
    diag("%s", float2str(123.f).c_str());
    diag("%s", float2str(1234.f).c_str());
    diag("%s", float2str(12345.f).c_str());
    diag("%s", float2str(123456.f).c_str());
    diag("%s", float2str(1234567.f).c_str());
    diag("%s", float2str(12345678.f).c_str());
    diag("%s", float2str(123456789.f).c_str());
    diag("time");
    diag("%s", time2str(5).c_str());
    diag("%s", time2str(35).c_str());
    diag("%s", time2str(59).c_str());
    diag("%s", time2str(60).c_str());
    diag("%s", time2str(61).c_str());
    diag("%s", time2str(601).c_str());
    diag("%s", time2str(3599).c_str());
    diag("%s", time2str(3600).c_str());
    diag("%s", time2str(3601).c_str());
    diag("%s", time2str(2123601).c_str());

    diag("Time: %s", local_time().c_str());
}

int main(int argc, char** argv)
{
    //plan_tests(15);
    plan_no_plan();

    test_proc1();

    return exit_status();
}
