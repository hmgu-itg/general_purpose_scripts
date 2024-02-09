#!/usr/bin/bash

cat <<EOF
<!DOCTYPE html>
<html>
<head>
  <title>JSON VIEW</title>
</head>
<style>
  #root {
      text-shadow: none;
      background: #D0D0D0;
      padding: 1em;
  }

  .renderjson a {
      text-decoration: none;
  }

  .renderjson .disclosure {
      color: crimson;
      font-size: 150%;
  }

  .renderjson .syntax {
      /* color: grey; */
      color: #D0D0D0;
  }

  .renderjson .string {
      color: red;
  }

  .renderjson .number {
      color: cyan;
  }

  .renderjson .boolean {
      color: plum;
  }

  .renderjson .key {
      color: darkblue;
  }

  .renderjson .keyword {
      color: lightgoldenrodyellow;
  }

  .renderjson .object.syntax {
      /* color: lightseagreen; */
      color: #D0D0D0;
  }

  .renderjson .array.syntax {
      color: lightsalmon;
  }
</style>

<body>
  <div id="root"></div>
  <script src="https://cdn.rawgit.com/caldwell/renderjson/master/renderjson.js"></script>

  <script>
    var data={
EOF

cat $1

cat <<EOF
}
    // Render toggable list in the container element
    document.getElementById("root").appendChild(
        renderjson(data)
    );
  </script> 
</body>
</html>

EOF
