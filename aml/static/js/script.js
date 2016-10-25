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
    $('#myTable').DataTable({
        "lengthMenu": [[10, 20, 30, 40, 50, -1], [10, 20, 30, 40, 50, "All"]],
        "scrollX":true
    });
    $('#myTable tbody').on('click', 'tr', function(){
        $(this).toggleClass('selected');
    });

});


