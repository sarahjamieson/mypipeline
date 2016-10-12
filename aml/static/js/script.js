$(document).ready(function() {
    /*
    $('.result').click(function() {
        $(this).hide();
        $('button').show();
    });
    $('button').hide();
    $('button').click(function(){
        $('.result').show();
    });
    */
    $('table').DataTable({
        "scrollX":true
    });
    $('button').click(function(){
        $('#overflow').removeClass("overflow");
    });
});
// button only works for first row.